"""Validation behavior: `check_global_coverage` partition checks (F5) and
`GriddedRegions` region-type validation (F6)."""

import numpy as np
import xarray as xr
import pytest

from regionate.integrate import check_global_coverage
from regionate import Region, GriddedRegions
from test_gridded_regions import initialize_spherical_grid


class _StubRegion:
    def __init__(self, name, mask):
        self.name = name
        self.mask = mask


def _regions(masks):
    class _Regions:
        pass
    obj = _Regions()
    obj.region_dict = {i: _StubRegion(str(i), m) for i, m in enumerate(masks)}
    return obj


def _blank(grid):
    return xr.zeros_like(grid._ds.geolon, dtype=bool)


def test_check_global_coverage_accepts_partition():
    grid = initialize_spherical_grid(N=6)
    west = _blank(grid); west.values[:, :3] = True
    east = ~west
    check_global_coverage(_regions([west, east]))   # exact partition -> no raise


def test_check_global_coverage_raises_on_gap():
    grid = initialize_spherical_grid(N=6)
    west = _blank(grid); west.values[:, :3] = True
    east = _blank(grid); east.values[:, 3:5] = True   # column 5 covered by nobody
    with pytest.raises(ValueError, match="no region"):
        check_global_coverage(_regions([west, east]))


def test_check_global_coverage_raises_on_overlap():
    grid = initialize_spherical_grid(N=6)
    west = _blank(grid); west.values[:, :4] = True
    east = _blank(grid); east.values[:, 3:] = True    # column 3 covered by both
    with pytest.raises(ValueError, match="more than one"):
        check_global_coverage(_regions([west, east]))


def test_griddedregions_rejects_bad_member_type():
    grid = initialize_spherical_grid(N=6)
    with pytest.raises(TypeError, match="Region"):
        GriddedRegions({"bad": object()}, grid)


def test_griddedregions_accepts_region():
    grid = initialize_spherical_grid(N=6)
    region = Region("r", np.array([0., 120., 240.]), np.array([0., 10., -10.]))
    gr = GriddedRegions({region.name: region}, grid)
    assert gr.region_dict["r"] is region
