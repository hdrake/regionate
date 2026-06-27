"""Load the ECCOv4r4 native lat-lon-cap (LLC90) grid as a symmetric xgcm.Grid
that regionate can consume.

The ECCO geometry file is distributed by NASA PO.DAAC and requires a (free)
NASA Earthdata Login. If the file is not already present under ``../data/`` it
is downloaded with ``earthaccess`` (which reads credentials from ``~/.netrc``;
run ``earthaccess.login(persist=True)`` once to set this up).

The native grid is MITgcm-staggered (vorticity points on the SW / 'left'
corner) and names its coordinates ``XC/YC`` (centers) and ``XG/YG`` (corners).
We rename those to the ``geolon*/geolat*`` convention and use
``sectionate.gridutils.symmetrize`` to build the equivalent symmetric ('outer')
grid, which regionate/sectionate operate on.
"""

import os
import numpy as np
import xarray as xr
import xgcm

ECCO_GEOMETRY_FILE = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
ECCO_GEOMETRY_SHORTNAME = "ECCO_L4_GEOMETRY_LLC0090GRID_V4R4"

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


def download_ECCO_geometry(data_dir="../data"):
    """Return the local path to the ECCO geometry file, downloading it from
    PO.DAAC via earthaccess if it is not already present."""
    path = os.path.join(data_dir, ECCO_GEOMETRY_FILE)
    if os.path.exists(path):
        return path
    import earthaccess
    earthaccess.login()  # uses ~/.netrc
    results = earthaccess.search_data(short_name=ECCO_GEOMETRY_SHORTNAME, count=1)
    earthaccess.download(results, local_path=data_dir)
    return path


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
    return xgcm.Grid(
        ds, periodic=False, autoparse_metadata=False,
        coords={"X": {"center": "i", "left": "i_g"},
                "Y": {"center": "j", "left": "j_g"}},
        face_connections=LLC90_FACE_CONNECTIONS,
    )


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
    """
    import regionmask
    ob = regionmask.defined_regions.natural_earth_v5_1_2.ocean_basins_50
    atl_ids = [int(n) for n, nm in zip(ob.numbers, ob.names) if nm in ATLANTIC_BASINS]
    facedim = grid._facedim
    lon, lat, depth = grid._ds["geolon"], grid._ds["geolat"], grid._ds["Depth"]
    arr = np.zeros(depth.shape, dtype=bool)
    for f in range(grid._ds.sizes[facedim]):
        ids = ob.mask(lon.isel({facedim: f}), lat.isel({facedim: f})).values
        arr[f] = np.isin(ids, atl_ids)
    return xr.DataArray(arr & (depth.values > 0), dims=depth.dims, coords=depth.coords)
