"""Boundary tracing on a single-tile grid that is periodic in X.

A region that wraps the periodic-X seam (geographically, the antimeridian) must trace
into a SINGLE stitched loop -- not split into two pieces at the seam -- and a strip that
circles the whole globe must trace into two latitude circles (each closed across the
seam). In every case the traced loops' velocity faces must reproduce the region's flux
convergence (the discrete divergence theorem)."""

import numpy as np
import xarray as xr
import xgcm

import sectionate as sec
from regionate.boundaries import grid_boundaries_from_mask


def periodic_grid(Nx=8, Ny=5):
    """A self-contained single-tile grid, periodic in X, MOM6 'symmetric' (outer)
    staggering. Longitude spans [-180, 180) so the periodic seam is the antimeridian."""
    xq = np.arange(Nx + 1); yq = np.arange(Ny + 1)
    xh = np.arange(Nx) + 0.5; yh = np.arange(Ny) + 0.5
    lonq = -180.0 + xq * (360.0 / Nx)
    lonh = -180.0 + xh * (360.0 / Nx)
    latq = -20.0 + yq * (40.0 / Ny)
    lath = -20.0 + yh * (40.0 / Ny)
    ds = xr.Dataset(coords={
        "xh": ("xh", xh), "xq": ("xq", xq.astype(float)),
        "yh": ("yh", yh), "yq": ("yq", yq.astype(float)),
        "geolon_c": (("yq", "xq"), np.broadcast_to(lonq[None, :], (Ny + 1, Nx + 1))),
        "geolat_c": (("yq", "xq"), np.broadcast_to(latq[:, None], (Ny + 1, Nx + 1))),
        "geolon": (("yh", "xh"), np.broadcast_to(lonh[None, :], (Ny, Nx))),
        "geolat": (("yh", "xh"), np.broadcast_to(lath[:, None], (Ny, Nx)))})
    return xgcm.Grid(
        ds, coords={"X": {"center": "xh", "outer": "xq"},
                    "Y": {"center": "yh", "outer": "yq"}},
        padding={"X": "periodic", "Y": "extend"}, autoparse_metadata=False)


def make_mask(grid, cells):
    arr = np.zeros_like(grid._ds.geolon.values, dtype=bool)
    for (j, i) in cells:
        arr[j, i] = True
    return xr.DataArray(arr, dims=grid._ds.geolon.dims, coords=grid._ds.geolon.coords)


def _convergence_and_boundary_flux(grid, mask, seed):
    """Cell-convergence vs boundary-flux for a random transport field whose periodic
    seam U-face is single-valued; the two must agree by the discrete divergence theorem
    (summing over every returned loop)."""
    Nx = grid._ds.sizes["xh"]; Ny = grid._ds.sizes["yh"]
    rng = np.random.default_rng(seed)
    U = rng.standard_normal((Ny, Nx + 1))          # (yh, xq)
    V = rng.standard_normal((Ny + 1, Nx))          # (yq, xh)
    U[:, Nx] = U[:, 0]                              # periodic seam U-face single-valued

    js, iss = np.where(mask.values)
    conv = 0.0
    for (jc, ic) in zip(js, iss):
        conv -= (U[jc, ic + 1] - U[jc, ic]) + (V[jc + 1, ic] - V[jc, ic])

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


def test_periodic_x_wrapping_band_stitches_into_one_loop():
    grid = periodic_grid(Nx=8, Ny=5)
    # a 2x2 block occupying the last and first columns -> adjacent across the seam
    mask = make_mask(grid, [(2, 7), (2, 0), (3, 7), (3, 0)])
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1                                  # single stitched loop
    assert np.array_equal(f_l[0], np.zeros_like(f_l[0]))                                 # single tile: one face, so zeros
    # the loop spans both sides of the seam (near lon -180 and +180)
    lons = np.asarray(lon_l[0])
    assert (lons < -90).any() and (lons > 90).any()


def test_periodic_x_wrapping_band_obeys_discrete_divergence_theorem():
    grid = periodic_grid(Nx=8, Ny=5)
    mask = make_mask(grid, [(2, 7), (2, 0), (3, 7), (3, 0)])
    conv, flux = _convergence_and_boundary_flux(grid, mask, seed=1)
    assert np.isclose(conv, flux, atol=1e-9)


def test_periodic_x_zonal_strip_is_two_latitude_circles():
    grid = periodic_grid(Nx=8, Ny=5)
    mask = make_mask(grid, [(2, i) for i in range(8)])    # full-width strip at j=2
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 2                                  # top + bottom latitude circle
    # neither circle carries a spurious radial seam segment: each is a constant-lat loop
    for k in range(2):
        assert np.allclose(np.asarray(lat_l[k]), np.asarray(lat_l[k])[0])
    conv, flux = _convergence_and_boundary_flux(grid, mask, seed=2)
    assert np.isclose(conv, flux, atol=1e-9)


def test_interior_region_away_from_seam_unaffected():
    grid = periodic_grid(Nx=8, Ny=5)
    mask = make_mask(grid, [(1, 3), (1, 4), (2, 3)])      # nowhere near the seam
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1
    conv, flux = _convergence_and_boundary_flux(grid, mask, seed=3)
    assert np.isclose(conv, flux, atol=1e-9)
