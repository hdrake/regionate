"""Load the ECCOv4r4 native lat-lon-cap (LLC90) grid as a symmetric xgcm.Grid
that regionate can consume, plus one month of native 3-D temperature fluxes.

The ECCO geometry and flux files are redistributed (with the ECCO Consortium's
permission) from the regionate/sectionate example archive on Zenodo,
https://doi.org/10.5281/zenodo.21051424 -- so no NASA Earthdata login is needed.
Files not already present under ``../data/`` are fetched on demand over HTTP,
mirroring ``load_example_model_grid.download_MOM6_example_data``.

The native grid is MITgcm-staggered (vorticity points on the SW / 'left'
corner) and names its coordinates ``XC/YC`` (centers) and ``XG/YG`` (corners).
We rename those to the ``geolon*/geolat*`` convention; sectionate consumes the
native 'left' staggering directly via ``corner_position``/``corner_offset``.
"""

import os
import urllib.request
import shutil
import numpy as np
import xarray as xr
import xgcm

# Concept DOI https://doi.org/10.5281/zenodo.21051424 resolves to this version
# record; files are served at <record>/files/<file_name>.
ZENODO_RECORD = "21051920"
ZENODO_FILES_URL = f"https://zenodo.org/records/{ZENODO_RECORD}/files/"

ECCO_GEOMETRY_FILE = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
# Monthly 3-D advective + diffusive temperature fluxes (2010-MM); MM in 1..12.
ECCO_TEMPERATURE_FLUX_FILE = (
    "OCEAN_3D_TEMPERATURE_FLUX_mon_mean_2010-{month:02d}_ECCO_V4r4_native_llc0090.nc"
)

# Canonical xgcm face_connections for the 13-tile LLC90 grid
# (from the xgcm ECCOv4 example), keyed by the ECCO 'tile' dimension.
LLC90_FACE_CONNECTIONS = {"tile": {
    0:  {"X": ((12, "Y", False), (3, "X", False)), "Y": (None,            (1, "Y", False))},
    1:  {"X": ((11, "Y", False), (4, "X", False)), "Y": ((0, "Y", False), (2, "Y", False))},
    2:  {"X": ((10, "Y", False), (5, "X", False)), "Y": ((1, "Y", False), (6, "X", False))},
    3:  {"X": ((0,  "X", False), (9, "Y", False)), "Y": (None,            (4, "Y", False))},
    4:  {"X": ((1,  "X", False), (8, "Y", False)), "Y": ((3, "Y", False), (5, "Y", False))},
    5:  {"X": ((2,  "X", False), (7, "Y", False)), "Y": ((4, "Y", False), (6, "Y", False))},
    6:  {"X": ((2,  "Y", False), (7, "X", False)), "Y": ((5, "Y", False), (10, "X", False))},
    7:  {"X": ((6,  "X", False), (8, "X", False)), "Y": ((5, "X", False), (10, "Y", False))},
    8:  {"X": ((7,  "X", False), (9, "X", False)), "Y": ((4, "X", False), (11, "Y", False))},
    9:  {"X": ((8,  "X", False), None),            "Y": ((3, "X", False), (12, "Y", False))},
    10: {"X": ((6,  "Y", False), (11, "X", False)), "Y": ((7, "Y", False), (2, "X", False))},
    11: {"X": ((10, "X", False), (12, "X", False)), "Y": ((8, "Y", False), (1, "X", False))},
    12: {"X": ((11, "X", False), None),            "Y": ((9, "Y", False), (0, "X", False))},
}}


def download_ECCO_data(file_name, data_dir="../data"):
    """Return the local path to an ECCO example file, downloading it from the
    Zenodo archive (https://doi.org/10.5281/zenodo.21051424) into ``data_dir``
    if it is not already present. No NASA Earthdata login is required."""
    destination_path = os.path.join(data_dir, file_name)
    if not os.path.exists(destination_path):
        os.makedirs(data_dir, exist_ok=True)
        print(f"File '{file_name}' being downloaded to {destination_path}.")
        # Stream to a temporary file and atomically rename on success, so an interrupted
        # transfer never leaves a truncated file that a later call mistakes for complete.
        tmp_path = destination_path + ".part"
        try:
            with urllib.request.urlopen(ZENODO_FILES_URL + file_name) as response, \
                    open(tmp_path, "wb") as out_file:
                shutil.copyfileobj(response, out_file)
            os.replace(tmp_path, destination_path)
        except BaseException:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            raise
        print(f"File '{file_name}' has completed download to {destination_path}.")
    else:
        print(f"File '{file_name}' already exists at {destination_path}. Skipping download.")
    return destination_path


def download_ECCO_geometry(data_dir="../data"):
    """Return the local path to the ECCO geometry file, downloading it from
    Zenodo if it is not already present."""
    return download_ECCO_data(ECCO_GEOMETRY_FILE, data_dir=data_dir)


# The LLC90 domain's southern edge folds onto itself, along the great circle
# through 65E and 115W under Antarctica. `face_connections` declares those four
# tile edges as `None` -- walls -- because its schema maps a whole tile edge to a
# whole tile edge, and this edge is glued to several partners over different index
# ranges: tile 0's runs into tile 3's, tile 9's and tile 12's in turn, and two of
# them fold onto themselves at the pivots. So one line of physical corners is
# stored twice with nothing in the metadata to say the two copies are the same
# point, and a region straddling it cannot close.
#
# Each block below is an affine map on the 'outer' corner lattice, `k -> const - k`
# along the named edge. Together they merge 179 corner pairs; every merged pair is
# bit-identically coincident in the archived coordinates (`validate_positions()`
# returns 0.0 m), which is a check on the declaration rather than its source.
LLC90_SOUTHERN_FOLD_BLOCKS = [
    #  (face, edge,     range,        face, edge,     const)
    (0,  "Y_low",  range(91),      3,  "Y_low",  116),
    (3,  "Y_low",  range(27),      3,  "Y_low",   26),
    (0,  "Y_low",  range(27),      9,  "X_high",  26),
    (9,  "X_high", range(91),      12, "X_high", 116),
    (12, "X_high", range(27),      12, "X_high",  26),
]


def llc90_southern_fold(n=91):
    """The LLC90 southern boundary fold, as explicit corner identifications."""
    from sectionate.topology import Identification

    def slot(f, edge, k):
        return (f, 0, k) if edge == "Y_low" else (f, k, n - 1)

    out = []
    for fa, ea, ks, fb, eb, const in LLC90_SOUTHERN_FOLD_BLOCKS:
        ks = [k for k in ks if 0 <= const - k < n]
        out.append(Identification(
            [slot(fa, ea, k) for k in ks],
            [slot(fb, eb, const - k) for k in ks],
            name=f"LLC90 southern fold {fa}.{ea} <-> {fb}.{eb}",
        ))
    return out


def load_ECCO_LLC90_grid(data_dir="../data"):
    """Load the ECCOv4r4 LLC90 grid as a native (MITgcm 'left'-staggered)
    ``xgcm.Grid``.

    The returned grid has tracer-center coordinates ``geolon``/``geolat`` and
    cell-corner coordinates ``geolon_c``/``geolat_c`` (on the native 'left'
    vorticity position; sectionate supports this directly via
    ``corner_position``/``corner_offset``), carries ``Depth`` for land/ocean
    masking, and encodes the LLC90 tile topology via ``face_connections`` -- so
    regionate can trace region boundaries across tile seams.
    """
    path = download_ECCO_geometry(data_dir=data_dir)
    ds = xr.open_dataset(path)
    ds = ds.rename({"XC": "geolon", "YC": "geolat",
                    "XG": "geolon_c", "YG": "geolat_c"})
    from sectionate.topology import declare_identifications

    grid = xgcm.Grid(
        ds, padding="fill", autoparse_metadata=False,
        coords={"X": {"center": "i", "left": "i_g"},
                "Y": {"center": "j", "left": "j_g"}},
        face_connections=LLC90_FACE_CONNECTIONS,
    )
    declare_identifications(grid, llc90_southern_fold())
    return grid


# Natural Earth ocean basins that make up the Atlantic sector (open Atlantic plus
# its embayments / marginal seas). The open ocean is split into several named
# basins -- notably the Sargasso Sea fills the centre of the North Atlantic -- so
# selecting only "*Atlantic*" would leave a large hole. The Mediterranean,
# Pacific and Arctic basins are deliberately excluded.
ATLANTIC_BASINS = (
    "North Atlantic Ocean", "South Atlantic Ocean", "Sargasso Sea",
    "Caribbean Sea", "Labrador Sea", "Gulf of Guinea", "Bay of Biscay",
    "Gulf of Saint Lawrence", "Gulf of Maine", "Bay of Fundy", "Ungava Bay",
    "Inner Seas", "Irish Sea", "Bristol Channel", "English Channel",
    "North Sea", "Río de la Plata", "Golfo San Jorge",
)


def atlantic_basin_mask(grid):
    """Boolean Atlantic-basin ocean mask on the ECCO grid, from the published
    Natural Earth ocean basins (``ATLANTIC_BASINS``) via ``regionmask``.

    Using published basin polygons (rather than a lon/lat box) keeps the mask
    geographically correct -- it follows the coastlines, fills the central
    Atlantic, and excludes the Pacific.

    The per-tile rasterization + stitching that a multi-tile grid needs is done by
    ``regionate.grid_conform.rasterize_per_tile`` (the same helper that lets
    ``regionate`` build boundary->mask regions on lat-lon-cap grids), rather than
    reimplemented here.
    """
    import regionmask
    from regionate.grid_conform import rasterize_per_tile
    ob = regionmask.defined_regions.natural_earth_v5_1_2.ocean_basins_50
    atl_ids = [int(n) for n, nm in zip(ob.numbers, ob.names) if nm in ATLANTIC_BASINS]

    def _atlantic(lon_2d, lat_2d):
        return np.isin(ob.mask(lon_2d, lat_2d).values, atl_ids)

    atlantic = rasterize_per_tile(grid, _atlantic)
    return atlantic & (grid._ds["Depth"] > 0)


def load_ECCO_temperature_flux(month=1, data_dir="../data"):
    """Load one month of ECCOv4r4 native 3-D temperature fluxes (default
    2010-01), downloading the file from Zenodo if needed.

    The returned single-time-step ``xarray.Dataset`` carries the advective heat
    fluxes ``ADVx_TH``/``ADVy_TH`` on the U/V velocity faces (units degC m3 s-1),
    indexed by the native vertical tracer dimension ``k`` (so a column integral
    is ``ADVx_TH.sum("k")``).

    The two flux components are promoted to float64. The archived data is
    float32, whose ~1e-7 relative precision (about 12 degC m3 s-1 on these
    O(1e8) transports) would otherwise dominate the discrete divergence-theorem
    residual and mask the exact boundary/interior closure.
    """
    path = download_ECCO_data(
        ECCO_TEMPERATURE_FLUX_FILE.format(month=month), data_dir=data_dir
    )
    ds = xr.open_dataset(path).isel(time=0)
    for var in ("ADVx_TH", "ADVy_TH"):
        ds[var] = ds[var].astype("float64")
    return ds
