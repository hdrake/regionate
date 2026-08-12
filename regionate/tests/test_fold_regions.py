"""Boundary tracing on a single-tile grid carrying a bipolar/tripolar north fold
(xgcm ``padding={..., "Y": {"fold": ...}}``, requires the north-fold boundary
released in xgcm >= 0.10.1).

The northern edge folds onto itself: top-row cell ``i`` is the fold-neighbour of
cell ``Nx-1-i`` (corner pivot, mirror about ``x=0``). A region straddling the fold
must trace into a SINGLE boundary loop with no spurious seam segment, and that
loop's velocity faces must reproduce the region's flux convergence exactly."""

import numpy as np
import xarray as xr
import xgcm
import pytest

import sectionate as sec
from regionate.boundaries import grid_boundaries_from_mask
from regionate import MaskRegions


def fold_grid(Nx=6, Ny=4):
    """A self-contained single-tile fold grid, MOM6 'symmetric' (outer) staggering,
    corner pivot.

    A cylinder, periodic in X, pinched shut as `j` rises: at row `j` the corners lie
    on an ellipse whose two semi-axes both shrink, and the top one has zero height,
    so the seam row collapses onto a segment traversed out and back. Corner `i` and
    its mirror are then exactly one point -- the corner-pivot fold identity -- with
    the two poles at `i = 0` and `i = Nx/2`.

    Both semi-axes shrink so the ellipses genuinely nest and every cell of the grid
    runs the same way round. A family whose axes moved in opposite directions would
    cross, leaving the cells nearest the seam inside out and their velocity faces
    signed backwards.
    """
    xq = np.arange(Nx + 1); yq = np.arange(Ny + 1)
    xh = np.arange(Nx) + 0.5; yh = np.arange(Ny) + 0.5

    def lonlat(j, i):
        theta = 2.0 * np.pi * np.asarray(i) / Nx
        t = np.asarray(j) / Ny
        return (25.0 * (1.0 - 0.5 * t) * np.cos(theta),
                60.0 + 18.0 * (1.0 - t) * np.sin(theta))

    LONc, LATc = lonlat(*np.meshgrid(yq, xq, indexing="ij"))
    LON, LAT = lonlat(*np.meshgrid(yh, xh, indexing="ij"))
    ds = xr.Dataset(coords={
        "xh": ("xh", xh), "xq": ("xq", xq.astype(float)),
        "yh": ("yh", yh), "yq": ("yq", yq.astype(float)),
        "geolon_c": (("yq", "xq"), LONc), "geolat_c": (("yq", "xq"), LATc),
        "geolon": (("yh", "xh"), LON), "geolat": (("yh", "xh"), LAT)})
    return xgcm.Grid(
        ds, coords={"X": {"center": "xh", "outer": "xq"},
                    "Y": {"center": "yh", "outer": "yq"}},
        padding={"X": "periodic", "Y": {"fold": "corner"}}, autoparse_metadata=False)


def make_mask(grid, cells):
    """cells: list of (j, i) center cells to set True."""
    arr = np.zeros_like(grid._ds.geolon.values, dtype=bool)
    for (j, i) in cells:
        arr[j, i] = True
    return xr.DataArray(arr, dims=grid._ds.geolon.dims, coords=grid._ds.geolon.coords)


def _convergence_and_boundary_flux(grid, mask, seed):
    """Return (cell-convergence, boundary-flux) for a random *fold-consistent* face
    transport field. The two must agree by the discrete divergence theorem."""
    Nx = grid._ds.sizes["xh"]; Ny = grid._ds.sizes["yh"]
    rng = np.random.default_rng(seed)
    U = rng.standard_normal((Ny, Nx + 1))          # (yh, xq)
    V = rng.standard_normal((Ny + 1, Nx))          # (yq, xh)
    # fold vector-sign constraint on the seam row: V[Ny, i] = -V[Ny, Nx-1-i]
    for i in range(Nx):
        mir = Nx - 1 - i
        if i < mir:
            a = 0.5 * (V[Ny, i] - V[Ny, mir]); V[Ny, i] = a; V[Ny, mir] = -a

    # convergence into the region = -sum of divergence over masked cells
    js, iss = np.where(mask.values)
    conv = 0.0
    for (jc, ic) in zip(js, iss):
        conv -= (U[jc, ic + 1] - U[jc, ic]) + (V[jc + 1, ic] - V[jc, ic])

    # boundary flux = sum of signed velocity faces of every traced loop
    # (single-tile sign convention: Usign=+1 if not Yinc; Vsign=+1 if Xinc)
    flux = 0.0
    i_l, j_l, f_l, _, _ = grid_boundaries_from_mask(grid, mask)
    for k in range(len(i_l)):
        uv = sec.uvindices_from_qindices(grid, i_l[k], j_l[k], f_c=f_l[k])
        for t in range(len(uv["var"])):
            if uv["var"][t] == "0":
                continue
            ii, jj = int(uv["i"][t]), int(uv["j"][t])
            if uv["var"][t] == "U":
                flux += (1 if not uv["Yinc"][t] else -1) * U[jj, ii]
            else:
                flux += (1 if uv["Xinc"][t] else -1) * V[jj, ii]
    return conv, flux


def test_fold_straddling_region_stitches_into_one_loop():
    grid = fold_grid(Nx=6, Ny=4)
    # a top-row cell and its across-fold mirror: connected ONLY through the fold
    mask = make_mask(grid, [(3, 1), (3, 4)])          # mirror of i=1 is Nx-1-1=4
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1                               # a single stitched loop
    # a single-tile grid is one face, so the face index is present and zero
    assert np.array_equal(f_l[0], np.zeros_like(f_l[0]))
    # its six boundary faces carry no seam face (the fold face is interior)
    uv = sec.uvindices_from_qindices(grid, i_l[0], j_l[0], f_c=f_l[0])
    assert sum(v != "0" for v in uv["var"]) == 6
    # and MaskRegions sees one region
    assert len(MaskRegions(mask, grid).region_dict) == 1


def test_fold_straddling_region_obeys_discrete_divergence_theorem():
    grid = fold_grid(Nx=6, Ny=4)
    mask = make_mask(grid, [(3, 1), (3, 4)])
    conv, flux = _convergence_and_boundary_flux(grid, mask, seed=1)
    assert np.isclose(conv, flux, atol=1e-9)


def test_single_fold_cell_keeps_seam_boundary_face():
    grid = fold_grid(Nx=6, Ny=4)
    # one top-row cell whose fold-mirror is OUTSIDE the mask: its north (fold) face
    # is a genuine boundary face and must be kept -> a 4-corner box, and the budget
    # still closes (sectionate signs the seam V-face correctly).
    mask = make_mask(grid, [(3, 1)])
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1
    uv = sec.uvindices_from_qindices(grid, i_l[0], j_l[0], f_c=f_l[0])
    assert sum(v != "0" for v in uv["var"]) == 4       # W, E, S, and the fold N face
    conv, flux = _convergence_and_boundary_flux(grid, mask, seed=2)
    assert np.isclose(conv, flux, atol=1e-9)


def test_fold_interior_region_unaffected():
    grid = fold_grid(Nx=6, Ny=4)
    # a region far from the seam traces as a plain box and closes the budget
    mask = make_mask(grid, [(1, 2), (1, 3), (2, 2)])
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1
    uv = sec.uvindices_from_qindices(grid, i_l[0], j_l[0], f_c=f_l[0])
    assert sum(v != "0" for v in uv["var"]) == 8       # L-tromino perimeter
    conv, flux = _convergence_and_boundary_flux(grid, mask, seed=3)
    assert np.isclose(conv, flux, atol=1e-9)
