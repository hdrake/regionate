"""Regression tests for issue #24: a 2D (X,Y) cell mask on a 3D (X,Y,Z) grid.

`MaskRegions(2d_mask, 3d_grid)` used to crash inside `xgcm.pad` with
`KeyError: "None of the DataArray's dims ('yh', 'xh') were found in axis coords."`
because `regionate.boundaries._pad_center` iterated over EVERY axis in `grid.axes`
(including the vertical Z axis) even though the center-point field it pads carries
only the horizontal tracer dims. A budget over a horizontal region is a perfectly
ordinary use of a 3D ocean grid, so the fix pads only the axes whose dims the field
actually has, and the seam-aware mask/boundary answer must be identical to the one on
the plain 2D grid (adding an unused Z axis changes nothing about horizontal topology).
"""

import numpy as np
import xarray as xr
import xgcm

import regionate
from regionate.boundaries import connected_components
from regionate.tests.test_gridded_regions import initialize_spherical_grid
from regionate.tests.test_connected_components import _mask_from_cells


def _add_z_axis(grid_2d):
    """Return a 3D `xgcm.Grid` identical to `grid_2d` but with an added vertical (Z)
    axis (`zl` centers, `zi` interfaces). The X/Y coords dict mirrors
    `initialize_spherical_grid` exactly, so the only difference is the extra Z axis."""
    zi = np.array([0.0, 10.0, 30.0, 60.0])       # interfaces (Z outer)
    zl = 0.5 * (zi[:-1] + zi[1:])                # centers (Z center)
    ds = grid_2d._ds.copy()
    ds = ds.assign_coords(
        zl=xr.DataArray(zl, dims=("zl",)),
        zi=xr.DataArray(zi, dims=("zi",)),
    )
    coords = {
        "X": {"outer": "xq", "center": "xh"},
        "Y": {"outer": "yq", "center": "yh"},
        "Z": {"center": "zl", "outer": "zi"},
    }
    return xgcm.Grid(
        ds, coords=coords,
        padding={"X": "periodic", "Y": "extend", "Z": "extend"},
        autoparse_metadata=False,
    )


def test_connected_components_ignores_unused_z_axis():
    """Lowest-level reproduction: `connected_components` on a 2D mask must give the same
    labeling on a 3D grid as on its 2D counterpart. Before the fix, the 3D call raised
    a KeyError inside `_pad_center`; the 2D call was always fine."""
    grid_2d = initialize_spherical_grid()
    grid_3d = _add_z_axis(grid_2d)

    # two clearly separated disks -> two components (cf. test_connected_components.py)
    mask = _mask_from_cells(grid_2d, [(1, 1), (4, 4)])

    labels_2d, ncomp_2d = connected_components(grid_2d, mask)
    labels_3d, ncomp_3d = connected_components(grid_3d, mask)  # used to raise KeyError

    assert ncomp_3d == ncomp_2d == 2
    assert np.array_equal(labels_3d.transpose(*labels_2d.dims).values, labels_2d.values)


def test_maskregions_ignores_unused_z_axis():
    """End-to-end: `MaskRegions(2d_mask, 3d_grid)` must succeed and produce the same
    region set as `MaskRegions(2d_mask, 2d_grid)`. Before the fix it raised a KeyError."""
    grid_2d = initialize_spherical_grid()
    grid_3d = _add_z_axis(grid_2d)

    mask = _mask_from_cells(grid_2d, [(1, 1), (4, 4)])

    regions_2d = regionate.MaskRegions(mask, grid_2d).region_dict
    regions_3d = regionate.MaskRegions(mask, grid_3d).region_dict  # used to raise

    assert len(regions_3d) == len(regions_2d) == 2
    for key in regions_2d:
        m2 = regions_2d[key].mask
        m3 = regions_3d[key].mask
        assert np.array_equal(m3.transpose(*m2.dims).values, m2.values)
