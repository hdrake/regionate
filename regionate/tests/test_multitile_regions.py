"""Boundary tracing and stitching across tile seams on multi-tile grids
(grids defined by xgcm `face_connections`, e.g. lat-lon-cap / cubed-sphere)."""

import numpy as np
import xarray as xr
import xgcm
import pytest

import sectionate as sec
from regionate.boundaries import grid_boundaries_from_mask
from regionate import MaskRegions


def rotated_two_tile_grid(Nc=4):
    """A native 'left'-staggered 2-tile grid whose tiles meet at a ROTATED seam:
    face 0's +X edge connects to face 1's Y axis (an X->Y, 90-degree connection).
    Each tile is given its own disjoint coordinate block (face 1 offset by +100),
    so the seam's two corner representations never coincide in lon/lat -- a region
    spanning the seam can therefore only be stitched via the grid topology
    (cell adjacency), not coordinate coincidence."""
    ng = Nc  # native 'left' corners are the same size as centers
    LONc = np.zeros((2, ng, ng)); LATc = np.zeros((2, ng, ng))
    LON = np.zeros((2, Nc, Nc));  LAT = np.zeros((2, Nc, Nc))
    for f in range(2):
        off = 100 * f
        LONc[f] = off + np.arange(ng)[None, :]
        LATc[f] = off + np.arange(ng)[:, None]
        LON[f] = off + np.arange(Nc)[None, :] + 0.5
        LAT[f] = off + np.arange(Nc)[:, None] + 0.5
    ds = xr.Dataset(coords={
        "i": ("i", np.arange(Nc)), "j": ("j", np.arange(Nc)),
        "i_g": ("i_g", np.arange(ng)), "j_g": ("j_g", np.arange(ng)),
        "face": ("face", [0, 1]),
        "geolon": (("face", "j", "i"), LON), "geolat": (("face", "j", "i"), LAT),
        "geolon_c": (("face", "j_g", "i_g"), LONc),
        "geolat_c": (("face", "j_g", "i_g"), LATc),
    })
    fc = {"face": {0: {"X": (None, (1, "Y", False))},
                   1: {"Y": ((0, "X", False), None)}}}
    grid_left = xgcm.Grid(
        ds, coords={"X": {"center": "i", "left": "i_g"},
                    "Y": {"center": "j", "left": "j_g"}},
        boundary="fill", fill_value=np.nan,
        face_connections=fc, autoparse_metadata=False,
    )
    return grid_left


def test_rotated_seam_region_stitches_into_one_loop():
    grid = rotated_two_tile_grid(Nc=4)
    # Both tiles fully in-mask: the rotated seam between them is internal, and
    # the region's outer boundary must be a single loop spanning both faces --
    # stitched across the rotated seam by topology alone (coords don't coincide).
    mask = xr.ones_like(grid._ds["geolon"]).astype(bool)
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1
    assert set(np.asarray(f_l[0]).tolist()) == {0, 1}


def two_face_grid(Nc=3):
    """Two faces side-by-side in longitude: face 0 spans [0, 90], face 1 [90, 180],
    joined at a single seam (face 0's right X edge -> face 1's left X edge). Carries
    both tracer-center (geolon/geolat) and corner (geolon_c/geolat_c) coordinates."""
    ng = Nc + 1
    yq = np.linspace(-45, 45, ng)
    yh = 0.5 * (yq[:-1] + yq[1:])
    lonq = [np.linspace(0, 90, ng), np.linspace(90, 180, ng)]
    lonh = [0.5 * (l[:-1] + l[1:]) for l in lonq]

    LONc = np.stack([np.broadcast_to(lonq[f], (ng, ng)) for f in range(2)])
    LATc = np.stack([np.broadcast_to(yq[:, None], (ng, ng)) for f in range(2)])
    LON = np.stack([np.broadcast_to(lonh[f], (Nc, Nc)) for f in range(2)])
    LAT = np.stack([np.broadcast_to(yh[:, None], (Nc, Nc)) for f in range(2)])

    ds = xr.Dataset(
        {},
        coords={
            "xq": (("xq",), np.arange(ng)), "yq": (("yq",), np.arange(ng)),
            "xh": (("xh",), np.arange(Nc)), "yh": (("yh",), np.arange(Nc)),
            "face": (("face",), [0, 1]),
            "geolon_c": (("face", "yq", "xq"), LONc),
            "geolat_c": (("face", "yq", "xq"), LATc),
            "geolon": (("face", "yh", "xh"), LON),
            "geolat": (("face", "yh", "xh"), LAT),
        },
    )
    fc = {"face": {0: {"X": (None, (1, "X", False))},
                   1: {"X": ((0, "X", False), None)}}}
    return xgcm.Grid(
        ds,
        coords={"X": {"outer": "xq", "center": "xh"},
                "Y": {"outer": "yq", "center": "yh"}},
        boundary="fill", fill_value=np.nan,
        face_connections=fc, autoparse_metadata=False,
    )


def make_mask(grid, cells):
    """cells: dict mapping face index -> list of (j, i) center cells to set True."""
    arr = np.zeros_like(grid._ds.geolon.values, dtype=bool)
    for f, lst in cells.items():
        for (j, i) in lst:
            arr[f, j, i] = True
    return xr.DataArray(arr, dims=grid._ds.geolon.dims, coords=grid._ds.geolon.coords)


def test_boundary_obeys_discrete_divergence_theorem():
    """A traced multi-tile boundary must obey the discrete divergence theorem:
    for any flux field, the net flux through the boundary's velocity faces equals
    the flux convergence summed over the masked cells. This is the property that
    makes regionate budgets consistent. Uses a region spanning the (non-rotated)
    tile seam and a seam-consistent synthetic transport field.

    (This synthetic grid is symmetric 'outer', so the shared seam U-face is stored
    on BOTH tiles and must be made single-valued by hand; on a native 'left' grid --
    e.g. real ECCO -- the seam face is stored once and any transport field works, as
    ``test_ecco_atlantic_basin_obeys_discrete_divergence_theorem`` checks.)"""
    grid = two_face_grid(Nc=6)
    Nc = grid._ds.sizes["xh"]; ng = Nc + 1
    # synthetic face transports; the shared seam U-face (face0 xq=Nc == face1 xq=0)
    # must be single-valued for the flux field to be physically consistent.
    umo = np.sin(np.arange(2 * Nc * ng).reshape(2, Nc, ng) * 0.07) + 0.3
    umo[1, :, 0] = umo[0, :, Nc]
    umo = xr.DataArray(umo, dims=("face", "yh", "xq"))
    vmo = xr.DataArray(np.cos(np.arange(2 * ng * Nc).reshape(2, ng, Nc) * 0.05) - 0.2,
                       dims=("face", "yq", "xh"))

    mask = make_mask(grid, {0: [(j, i) for j in range(1, 5) for i in range(3, Nc)],
                            1: [(j, i) for j in range(1, 5) for i in range(0, 3)]})
    convergence = float((-(grid.diff(umo, "X") + grid.diff(vmo, "Y"))).where(mask, 0.).sum())

    U, V = umo.values, vmo.values
    flux = 0.0
    for region in MaskRegions(mask, grid).region_dict.values():
        uv = sec.uvindices_from_qindices(grid, region.i_c, region.j_c, f_c=region.f_c)
        for k in range(len(uv["var"])):
            if uv["var"][k] == "0":
                continue
            f, i, j = int(uv["face"][k]), int(uv["i"][k]), int(uv["j"][k])
            flux += int(uv["Lsign"][k]) * (U[f, j, i] if uv["var"][k] == "U" else V[f, j, i])

    assert np.isclose(convergence, flux, atol=1e-9)


def test_seam_spanning_region_stitches_into_one_loop():
    grid = two_face_grid(Nc=3)
    # East column of face 0 + west column of face 1: a strip straddling the seam.
    mask = make_mask(grid, {0: [(0, 2), (1, 2), (2, 2)],
                            1: [(0, 0), (1, 0), (2, 0)]})
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)

    assert len(i_l) == 1                                  # a single stitched loop
    faces = set(np.asarray(f_l[0]).tolist())
    assert faces == {0, 1}                                # spans both faces
    # The seam faces are internal: the duplicated seam corners carry no velocity
    # face, so sectionate returns fewer faces than there are corners.
    lon_uv, lat_uv = sec.uvcoords_from_qindices(grid, i_l[0], j_l[0], f_c=f_l[0])
    assert len(lon_uv) == 10
    assert len(lon_l[0]) == 12


def test_seam_terminating_region_keeps_seam_edge():
    grid = two_face_grid(Nc=3)
    # East column of face 0 only: its east neighbour (face 1) is outside the mask,
    # so the seam edge at lon=90 is a real boundary and must be kept.
    mask = make_mask(grid, {0: [(0, 2), (1, 2), (2, 2)]})
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)

    assert len(i_l) == 1
    assert set(np.asarray(f_l[0]).tolist()) == {0}       # never leaves face 0
    assert np.isclose(np.asarray(lon_l[0]), 90.).any()   # seam edge retained
    # No internal seam faces here: every corner contributes a velocity face.
    lon_uv, _ = sec.uvcoords_from_qindices(grid, i_l[0], j_l[0], f_c=f_l[0])
    assert len(lon_uv) == len(lon_l[0])


def test_interior_cell_single_face_box():
    grid = two_face_grid(Nc=3)
    mask = make_mask(grid, {1: [(1, 1)]})
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1
    assert set(np.asarray(f_l[0]).tolist()) == {1}
    assert len(lon_l[0]) == 4                             # a 4-corner cell box


def test_maskregions_threads_face_index():
    grid = two_face_grid(Nc=3)
    mask = make_mask(grid, {0: [(0, 2), (1, 2), (2, 2)],
                            1: [(0, 0), (1, 0), (2, 0)]})
    regions = MaskRegions(mask, grid).region_dict
    assert len(regions) == 1
    region = regions[0]
    assert region.f_c is not None
    assert set(np.asarray(region.f_c).tolist()) == {0, 1}
