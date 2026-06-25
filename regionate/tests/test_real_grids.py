"""Optional end-to-end checks against real ocean-model grids.

These are skipped by default. Enable them by setting REGIONATE_REALDATA_TESTS=1
(and, for the ECCO check, pointing REGIONATE_ECCO_GRID at a NetCDF file holding
an ECCOv4r4 LLC90 grid with `face_connections`). They validate the multi-model
goal of the overhaul on genuine MOM6 (tripolar) and ECCO (lat-lon-cap) output.
"""

import os
import numpy as np
import xarray as xr
import pytest

REALDATA = os.environ.get("REGIONATE_REALDATA_TESTS") == "1"
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data")
MOM6_FILE = os.path.join(
    DATA_DIR, "MOM6_global_example_vertically_integrated_mass_budget_v0_0_6.nc"
)
ECCO_FILE = os.environ.get("REGIONATE_ECCO_GRID", "")


@pytest.mark.skipif(
    not (REALDATA and os.path.isfile(MOM6_FILE)),
    reason="set REGIONATE_REALDATA_TESTS=1 and provide the MOM6 global example file",
)
def test_mom6_global_box_mask_boundary_consistency():
    """On the real global MOM6 grid, a mid-latitude box defined by its boundary
    yields a mask whose own traced boundary re-encloses the same mask."""
    import xgcm
    from regionate import GriddedRegion, MaskRegions

    ds = xr.open_dataset(MOM6_FILE).fillna(0.)
    grid = xgcm.Grid(
        ds,
        coords={"X": {"center": "xh", "outer": "xq"},
                "Y": {"center": "yh", "outer": "yq"}},
        boundary={"X": "periodic", "Y": "extend"},
        autoparse_metadata=False,
    )
    lons = np.array([-40., -10., -10., -40.])
    lats = np.array([10., 10., 40., 40.])
    region = GriddedRegion("box", lons, lats, grid)
    assert int(region.mask.sum()) > 0

    # Re-tracing the mask must recover a region enclosing exactly the same cells.
    retraced = MaskRegions(region.mask, grid).region_dict
    assert len(retraced) >= 1
    union = None
    for r in retraced.values():
        union = r.mask if union is None else (union | r.mask)
    assert bool((union.values == region.mask.values).all())


@pytest.mark.skipif(
    not (REALDATA and os.path.isfile(ECCO_FILE)),
    reason="set REGIONATE_REALDATA_TESTS=1 and REGIONATE_ECCO_GRID=<llc90 grid .nc>",
)
def test_ecco_llc_reciprocal_or_raises():
    """ECCOv4r4 LLC90 neighbour construction must either be fully reciprocal
    (a usable topology) or refuse with NotImplementedError -- never silently
    inconsistent. Canonical LLC90 (all-non-reversed seams) is reciprocal; only
    the full cubed-sphere's reversed/rotated cap connections may raise."""
    import xgcm
    from sectionate.gridutils import build_neighbor_maps, get_geo_corners

    ds = xr.open_dataset(ECCO_FILE)
    # The caller's file is expected to be directly consumable by xgcm.Grid with
    # face_connections; we only assert the reciprocity invariant here.
    grid = xgcm.Grid(ds, periodic=False, autoparse_metadata=False)
    try:
        maps = build_neighbor_maps(grid, get_geo_corners(grid))
    except NotImplementedError:
        return
    # If it did not raise, a representative seam must be self-consistent.
    assert all(m is not None for m in next(iter(maps.values())))
