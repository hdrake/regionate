"""The grid topology used to trace a region's boundary can be chosen by the caller.

Up to v0.5.7 every boundary was traced with ``topology="MOM-tripolar"``, whatever the grid.
That is still the default, but a grid without a tripolar fold can now say so.
"""
import numpy as np
import pytest
import xarray as xr
import xgcm

import sectionate
from regionate import GriddedRegion
from regionate.tests.test_gridded_regions import initialize_spherical_grid

LONS = np.array([0., 120., 240., 360.])
LATS = np.array([0., 0., 0., 0.])
BOX_LONS = np.array([100., 160., 160., 100.])      # on the tripolar test grid, south of its cap
BOX_LATS = np.array([40., 40., 50., 50.])

# sectionate reads how a tripolar grid folds from its coordinates from this version on;
# earlier versions pair the points across the fold incorrectly (MOM6-community/sectionate#53)
fold_aware = hasattr(sectionate.section, "fold_partner_offset")


@pytest.fixture
def recorded_topologies(monkeypatch):
    seen = []
    grid_section = sectionate.grid_section

    def recording(grid, lons, lats, **kwargs):
        seen.append(kwargs.get("topology"))
        return grid_section(grid, lons, lats, **kwargs)

    monkeypatch.setattr(sectionate, "grid_section", recording)
    return seen


def test_default_topology_is_tripolar(recorded_topologies):
    grid = tripolar_grid()
    region = GriddedRegion("r", BOX_LONS, BOX_LATS, grid)
    assert region.topology == "MOM-tripolar"
    assert recorded_topologies == ["MOM-tripolar"]


def test_latlon_topology_is_passed_to_sectionate(recorded_topologies):
    region = GriddedRegion("r", LONS, LATS, initialize_spherical_grid(), topology="latlon")
    assert region.topology == "latlon"
    assert recorded_topologies == ["latlon"]


def test_topology_makes_no_difference_away_from_the_fold():
    grid = tripolar_grid()
    region = GriddedRegion("r", BOX_LONS, BOX_LATS, grid, topology="MOM-tripolar")
    same = GriddedRegion("r", BOX_LONS, BOX_LATS, grid, topology="latlon")
    assert np.array_equal(region.i_c, same.i_c) and np.array_equal(region.j_c, same.j_c)
    assert np.array_equal(region.mask, same.mask)


@pytest.mark.skipif(not fold_aware, reason="needs a sectionate that fixes MOM6-community/sectionate#53")
def test_default_topology_on_a_grid_without_a_fold_is_an_error():
    with pytest.raises(ValueError, match="does not fold onto itself"):
        GriddedRegion("r", LONS, LATS, initialize_spherical_grid())
    GriddedRegion("r", LONS, LATS, initialize_spherical_grid(), topology="latlon")


@pytest.mark.skipif(not fold_aware, reason="needs a sectionate that fixes MOM6-community/sectionate#53")
def test_boundary_across_the_fold():
    grid = tripolar_grid()
    lons = np.arange(0., 360., 30.) + 3.
    lats = np.full(lons.size, 75.)                  # a latitude circle inside the bipolar cap
    mask = xr.ones_like(grid._ds.geolon, dtype=bool)   # not under test; skips the polygon mask
    region = GriddedRegion("cap", lons, lats, grid, mask=mask, topology="MOM-tripolar")

    # every step joins adjacent corners, or the two indices of one point on the fold
    lon_c, lat_c = grid._ds.geolon_c.values, grid._ds.geolat_c.values
    i_c, j_c = region.i_c, region.j_c
    assert j_c.max() == lon_c.shape[0] - 1          # the boundary reaches the fold
    di = np.abs(np.diff(i_c)); di = np.minimum(di, NX - di)
    dj = np.abs(np.diff(j_c))
    twins = (di + dj) != 1
    assert twins.any()
    d = sectionate.distance_on_unit_sphere(
        lon_c[j_c[:-1], i_c[:-1]], lat_c[j_c[:-1], i_c[:-1]],
        lon_c[j_c[1:], i_c[1:]], lat_c[j_c[1:], i_c[1:]])
    assert np.all(d[twins] < 1.)                    # metres

    # one velocity face per step between distinct corners, so fewer than corners - 1
    assert len(region.lons_uv) == i_c.size - 1 - twins.sum()


# ---- an analytic MOM-style tripolar grid, adapted from MOM6-community/sectionate#54 ----
NX, NYC, LAT_SOUTH, LAT0, DLAT, LON0 = 72, 12, 30., 60., 2., 80.


def tripolar_grid():
    """Regular in latitude and longitude south of LAT0; north of it a bipolar cap whose poles sit
    on vorticity points 0 and NX/2 of its southern edge. The top row of vorticity points is the
    fold: (top, i) and (top, NX - i) are the same point."""
    m = 2*NX
    lam = 2*np.pi*np.arange(m+1)/m
    lat_s = np.arange(LAT_SOUTH, LAT0 + DLAT/4, DLAT/2)
    x_s = np.broadcast_to(LON0 + np.rad2deg(lam), (lat_s.size, m+1))
    y_s = np.broadcast_to(lat_s[:, None], (lat_s.size, m+1))
    a = np.tan(np.deg2rad(90. - LAT0)/2)
    with np.errstate(divide="ignore"):
        tau = np.log(np.abs(1/np.tan(lam/2)))
    frac = np.arange(1, 2*NYC+1)[:, None]/(2*NYC)
    sigma = np.where(lam < np.pi, -1., 1.)[None, :]*(np.pi/2)*(1. + frac)
    with np.errstate(over="ignore", invalid="ignore"):
        e = np.exp(tau[None, :] + 1j*sigma)
        z = a*(e + 1)/(e - 1)
    z[-1] = a*np.tanh(tau/2)
    poles = np.isinf(tau)
    z[:, poles & (tau > 0)] = a
    z[:, poles & (tau < 0)] = -a
    x = np.concatenate([x_s, LON0 + np.rad2deg(np.angle(z))], axis=0)
    y = np.concatenate([y_s, 90. - 2*np.rad2deg(np.arctan(np.abs(z)))], axis=0)
    ny = (x.shape[0] - 1)//2
    ds = xr.Dataset(coords={
        "xh": np.arange(NX) + 0.5, "xq": np.arange(NX+1) + 0.,
        "yh": np.arange(ny) + 0.5, "yq": np.arange(ny+1) + 0.,
        "geolon_c": (("yq", "xq"), x[::2, ::2]), "geolat_c": (("yq", "xq"), y[::2, ::2]),
        "geolon": (("yh", "xh"), x[1::2, 1::2]), "geolat": (("yh", "xh"), y[1::2, 1::2]),
    })
    return xgcm.Grid(
        ds, coords={"X": {"center": "xh", "outer": "xq"}, "Y": {"center": "yh", "outer": "yq"}},
        boundary={"X": "periodic", "Y": "extend"}, autoparse_metadata=False,
    )
