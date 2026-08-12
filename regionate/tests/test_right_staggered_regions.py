"""
Boundaries on a single-tile grid whose vorticity sits at the 'right' position.

Every other single-tile fixture in this suite is 'outer'-staggered, so nothing
exercised the corner frame the tracer works in. It used to depend on the
staggering: a 'right' grid was traced one corner off, which split a region
wrapping the periodic seam into two loops with a spurious seam face and a budget
that did not close. The corner topology is indexed on the 'outer' lattice whatever
the native staggering, so the frame is now the same for every grid.
"""

import numpy as np
import pytest
import xarray as xr
import xgcm

import sectionate as sec
from regionate.boundaries import grid_boundaries_from_mask


def right_grid(Nx=6, Ny=4):
    """A single-tile 'right'-staggered grid, periodic in X and walled in Y."""
    xq = (np.arange(Nx) + 1) * (360.0 / Nx)
    yq = np.linspace(-40.0, 40.0, Ny + 1)[1:]
    xh = xq - (360.0 / Nx) / 2.0
    yh = yq - (80.0 / Ny) / 2.0
    lon_c, lat_c = np.meshgrid(xq, yq)
    lon, lat = np.meshgrid(xh, yh)
    ds = xr.Dataset(coords={
        "xh": ("xh", xh), "yh": ("yh", yh), "xq": ("xq", xq), "yq": ("yq", yq),
        "geolon_c": (("yq", "xq"), lon_c), "geolat_c": (("yq", "xq"), lat_c),
        "geolon": (("yh", "xh"), lon), "geolat": (("yh", "xh"), lat)})
    return xgcm.Grid(
        ds, coords={"X": {"center": "xh", "right": "xq"},
                    "Y": {"center": "yh", "right": "yq"}},
        padding={"X": "periodic", "Y": "extend"}, autoparse_metadata=False)


def _mask(grid, cells, Nx=6, Ny=4):
    m = np.zeros((Ny, Nx), dtype=bool)
    for j, i in cells:
        m[j, i] = True
    return xr.DataArray(m, dims=("yh", "xh"))


def _budget_closes(grid, mask, i_c, j_c, f_c):
    """The discrete divergence theorem: the flux through the loop equals the
    convergence into the cells it encloses."""
    rng = np.random.default_rng(0)
    nyc, nxc = mask.shape
    u = rng.normal(size=(nyc, nxc))          # 'right': velocities are (Nc, Nc)
    v = rng.normal(size=(nyc, nxc))
    ds = grid._ds.assign({"ut": (("yh", "xq"), u), "vt": (("yq", "xh"), v)})
    g = xgcm.Grid(ds, coords={a: dict(grid.axes[a].coords) for a in ("X", "Y")},
                  padding={a: grid.axes[a].padding for a in ("X", "Y")},
                  autoparse_metadata=False)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = sec.convergent_transport(g, i_c, j_c, f_c=f_c, utr="ut", vtr="vt",
                                       positive_in=mask)
    flux = float(out["conv_mass_transport"].sum().values)
    # convergence straight from the staggered arrays, over the masked cells
    conv = 0.0
    for j in range(nyc):
        for i in range(nxc):
            if not bool(mask.values[j, i]):
                continue
            conv += u[j, i - 1] - u[j, i]            # X periodic: i-1 wraps
            conv += (v[j - 1, i] if j > 0 else 0.0) - v[j, i]
    return flux, conv


def test_interior_region_closes_its_budget():
    grid = right_grid()
    mask = _mask(grid, [(1, 2), (1, 3), (2, 2), (2, 3)])
    i_l, j_l, f_l, _, _ = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1
    flux, conv = _budget_closes(grid, mask, i_l[0], j_l[0], f_l[0])
    assert flux == pytest.approx(conv, abs=1e-12)


def test_a_region_wrapping_the_periodic_seam_is_one_loop_and_closes():
    """
    The case the corner frame used to get wrong. A band spanning the seam is one
    region with one boundary, and the seam is interior to it -- not an edge of it.
    """
    grid = right_grid()
    mask = _mask(grid, [(1, 5), (1, 0), (2, 5), (2, 0)])
    i_l, j_l, f_l, _, _ = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1, "a seam-wrapping region should trace as ONE loop"

    uv = sec.uvindices_from_qindices(grid, i_l[0], j_l[0], f_c=f_l[0])
    faces = list(zip(uv["var"].tolist(), uv["i"].tolist(), uv["j"].tolist()))
    assert len(faces) == len(set(faces)) == 8, "8 faces bound a 2x2 block"

    flux, conv = _budget_closes(grid, mask, i_l[0], j_l[0], f_l[0])
    assert flux == pytest.approx(conv, abs=1e-12)
