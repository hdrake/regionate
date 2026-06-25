"""Boundary tracing and stitching across tile seams on multi-tile grids
(grids defined by xgcm `face_connections`, e.g. lat-lon-cap / cubed-sphere)."""

import numpy as np
import xarray as xr
import xgcm
import pytest

import sectionate as sec
from regionate.boundaries import grid_boundaries_from_mask
from regionate import MaskRegions


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
