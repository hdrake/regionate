"""Topology-aware connected-component labeling of a cell mask.

`regionate.boundaries.connected_components` must label two in-mask cells into the
same component exactly when the grid's own topology makes them edge-neighbours --
across periodic axes, the north fold, and multi-tile `face_connections` -- and
must NOT merge cells that are only planar-adjacent across a wall. Each test pairs
the seam-aware answer against the (wrong) answer a planar labeller would give.
"""

import numpy as np
import xarray as xr

from regionate.boundaries import connected_components
from regionate.tests.test_gridded_regions import initialize_spherical_grid
from regionate.tests.test_multitile_regions import two_face_grid, make_mask as make_mask_mt
from regionate.tests.test_fold_regions import fold_grid, make_mask as make_mask_fold


def _mask_from_cells(grid, cells):
    arr = np.zeros_like(grid._ds.geolon.values, dtype=bool)
    for (j, i) in cells:
        arr[j, i] = True
    return xr.DataArray(arr, dims=grid._ds.geolon.dims, coords=grid._ds.geolon.coords)


def test_labels_align_with_mask_and_empty():
    grid = initialize_spherical_grid()
    empty = _mask_from_cells(grid, [])
    labels, ncomp = connected_components(grid, empty)
    assert ncomp == 0
    assert bool((labels == -1).all())

    mask = _mask_from_cells(grid, [(1, 1), (1, 2)])
    labels, ncomp = connected_components(grid, mask)
    # labels are >= 0 exactly on the mask, -1 elsewhere
    assert bool(((labels >= 0) == mask).all())
    assert set(np.unique(labels.values[mask.values]).tolist()) == set(range(ncomp))


def test_separated_disks_are_distinct_components():
    grid = initialize_spherical_grid()
    mask = _mask_from_cells(grid, [(1, 1), (4, 4)])  # far apart, no shared edge
    _, ncomp = connected_components(grid, mask)
    assert ncomp == 2


def test_periodic_seam_merges_wrapped_cells():
    """Cells at i=0 and i=N-1 in the same row touch ONLY across the periodic-X seam:
    seam-aware labeling merges them (1 component); a planar labeller would see 2."""
    grid = initialize_spherical_grid()
    N = grid._ds.sizes["xh"]
    mask = _mask_from_cells(grid, [(2, 0), (2, N - 1)])
    _, ncomp = connected_components(grid, mask)
    assert ncomp == 1


def test_global_strip_is_one_component_and_complement_is_two():
    """The zonal band |lat|<=10 circles the globe -> a single seam-wrapping component.
    Its complement is a disconnected north band and south band -> two components. This
    is exactly the region/loop distinction MaskRegions is being rebuilt around."""
    grid = initialize_spherical_grid()
    yh = grid._ds.yh
    strip = (np.abs(yh) <= 10).broadcast_like(grid._ds.geolon).astype(bool)
    _, ncomp = connected_components(grid, strip)
    assert ncomp == 1

    _, ncomp_inv = connected_components(grid, ~strip)
    assert ncomp_inv == 2


def test_fold_seam_merges_mirror_cells():
    """The two top-row cells (3,1) and (3,4)=(3,Nx-1-1) are mirror-neighbours across
    the bipolar fold; seam-aware labeling merges them, a planar labeller would not."""
    grid = fold_grid(Nx=6, Ny=4)
    mask = make_mask_fold(grid, [(3, 1), (3, 4)])
    _, ncomp = connected_components(grid, mask)
    assert ncomp == 1

    # a top-row cell whose mirror is NOT in the mask stays its own component
    solo = make_mask_fold(grid, [(3, 1), (0, 1)])  # (0,1) is interior, not adjacent
    _, ncomp_solo = connected_components(grid, solo)
    assert ncomp_solo == 2


def test_multitile_seam_merges_across_faces():
    """Face 0's east column and face 1's west column meet only across the tile seam:
    seam-aware labeling merges them into one component spanning both faces."""
    grid = two_face_grid(Nc=3)
    seam = make_mask_mt(grid, {0: [(0, 2), (1, 2), (2, 2)],
                               1: [(0, 0), (1, 0), (2, 0)]})
    _, ncomp = connected_components(grid, seam)
    assert ncomp == 1

    # cells buried in each face's interior, not sharing the seam -> two components
    apart = make_mask_mt(grid, {0: [(1, 0)], 1: [(1, 2)]})
    _, ncomp_apart = connected_components(grid, apart)
    assert ncomp_apart == 2
