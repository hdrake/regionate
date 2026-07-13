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


def test_overlap_alignment_keeps_indices_consistent_with_coords():
    """F3: aligning region boundaries to their shared overlap section must roll the
    grid indices (i_c/j_c/f_c) in lockstep with the coordinates, so each corner's
    stored index still points at that corner's coordinate (they were left stale)."""
    from regionate import GriddedRegion, GriddedRegions
    from regionate.overlaps import align_boundaries_with_overlap_sections

    grid = initialize_spherical_grid(N=12)
    gc = grid._ds.geolon_c.values
    gcl = grid._ds.geolat_c.values
    # two adjacent boxes sharing the lon=120 edge (names must be int-parseable)
    A = GriddedRegion("0", np.array([60., 120., 120., 60.]), np.array([-30., -30., 30., 30.]), grid)
    B = GriddedRegion("1", np.array([120., 180., 180., 120.]), np.array([-30., -30., 30., 30.]), grid)
    regs = GriddedRegions({"0": A, "1": B}, grid)
    regs.find_all_overlaps()
    assert ("0", "1") in regs.overlaps                 # they share the lon=120 edge

    def indices_match_coords(r):
        n = len(r.lons_c)
        lon = gc[np.asarray(r.j_c), np.asarray(r.i_c)][:n]
        lat = gcl[np.asarray(r.j_c), np.asarray(r.i_c)][:n]
        return (np.allclose(np.mod(lon, 360.), np.mod(r.lons_c, 360.))
                and np.allclose(lat, r.lats_c))

    align_boundaries_with_overlap_sections(regs)
    for name in ("0", "1"):
        assert indices_match_coords(regs.region_dict[name])


def test_gr_child_roundtrip_preserves_gridded_indices(tmp_path):
    """F7: `.gr` children must round-trip as gridded sections carrying their stored
    i_c/j_c/f_c -- open_gr previously rebuilt them as bare `sec.Section` from coords
    only, discarding the indices."""
    import xgcm
    import sectionate as sec
    from regionate import GriddedRegion
    from regionate.region import open_gr

    grid = initialize_spherical_grid(N=12)
    region = GriddedRegion("reg", np.array([60., 120., 120., 60.]),
                           np.array([-30., -30., 30., 30.]), grid)
    child = sec.GriddedSection(
        sec.Section("child", sec.coords_from_lonlat(region.lons_c[:4], region.lats_c[:4])),
        grid, i_c=region.i_c[:4], j_c=region.j_c[:4], f_c=None)
    region.children = {"child": child}
    region.to_gr(str(tmp_path))

    def ds_to_grid(ds):
        return xgcm.Grid(ds, coords={'X': {'outer': 'xq', 'center': 'xh'},
                                     'Y': {'outer': 'yq', 'center': 'yh'}},
                         padding={"X": "periodic", "Y": "extend"}, autoparse_metadata=False)

    reloaded = open_gr(f"{tmp_path}/reg.gr", ds_to_grid)
    cr = reloaded.children["child"]
    assert isinstance(cr, sec.GriddedSection)
    assert np.array_equal(np.asarray(cr.i_c), np.asarray(region.i_c[:4]))
    assert np.array_equal(np.asarray(cr.j_c), np.asarray(region.j_c[:4]))
