import numpy as np
from shapely.geometry import Polygon, MultiPolygon

from regionate.geometry import (
    normalize_lon,
    has_antimeridian_crossing,
    polygon_to_lon360,
    split_at_antimeridian,
)
from regionate.grid_conform import mask_from_grid_boundaries

# Reuse the single-tile spherical-grid fixture from the gridded-regions tests.
from regionate.tests.test_gridded_regions import initialize_spherical_grid


# ---------------------------------------------------------------------------
# geometry.py unit tests
# ---------------------------------------------------------------------------

def _all_lons_within(geom, lo=-180.0, hi=180.0, tol=1e-9):
    polys = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
    for p in polys:
        lons = np.array([c[0] for c in p.exterior.coords])
        assert np.all(lons >= lo - tol) and np.all(lons <= hi + tol)


def test_normalize_lon_roundtrip():
    lons = np.array([-300., -180., -1., 0., 1., 180., 359., 360., 540.])
    once = normalize_lon(lons)
    twice = normalize_lon(once)
    assert np.all(once >= -180.) and np.all(once < 180.)
    assert np.allclose(once, twice)
    # Equivalent longitudes map to the same normalized value.
    assert np.isclose(normalize_lon(370.), normalize_lon(10.))
    assert np.isclose(normalize_lon(-170.), normalize_lon(190.))


def test_has_antimeridian_crossing():
    crossing = Polygon([(170, -10), (-170, -10), (-170, 10), (170, 10)])
    interior = Polygon([(10, -10), (40, -10), (40, 10), (10, 10)])
    assert has_antimeridian_crossing(crossing)
    assert not has_antimeridian_crossing(interior)


def test_polygon_to_lon360():
    crossing = Polygon([(170, -10), (-170, -10), (-170, 10), (170, 10)])
    g360 = polygon_to_lon360(crossing)
    lons = np.array([c[0] for c in g360.exterior.coords])
    assert np.all(lons >= 0.) and np.all(lons <= 360.)
    # In [0,360] space the ring no longer has a huge discontinuity.
    assert not has_antimeridian_crossing(g360)


def test_split_dateline_crossing_into_multipolygon():
    crossing = Polygon([(170, -10), (-170, -10), (-170, 10), (170, 10)])
    split = split_at_antimeridian(crossing)
    assert split.geom_type == "MultiPolygon"
    assert len(split.geoms) == 2
    _all_lons_within(split)
    # The two pieces should hug the dateline (one near +180, one near -180).
    bounds = sorted(p.bounds for p in split.geoms)
    assert any(b[2] >= 179.0 for b in bounds)   # a piece reaching +180
    assert any(b[0] <= -179.0 for b in bounds)  # a piece reaching -180


def test_split_noncrossing_passthrough():
    interior = Polygon([(10, -10), (40, -10), (40, 10), (10, 10)])
    split = split_at_antimeridian(interior)
    assert split.geom_type == "Polygon"
    _all_lons_within(split)
    # Area is preserved for a simple non-crossing box.
    assert np.isclose(split.area, interior.area)


# ---------------------------------------------------------------------------
# Functional tests against the spherical grid fixture
# ---------------------------------------------------------------------------

def test_mask_interior_box():
    grid = initialize_spherical_grid(N=6)
    # Grid centers sit at xh = 30,90,...,330 ; yh = -50,-30,...,50.
    # A box [60,120] x [-20,20] should select exactly the xh=90, yh in {-10,10} cells.
    lons = np.array([60., 120., 120., 60., 60.])
    lats = np.array([-20., -20., 20., 20., -20.])
    mask = mask_from_grid_boundaries(lons, lats, grid)

    geolon = grid._ds["geolon"]
    geolat = grid._ds["geolat"]
    expected = (geolon == 90.) & (np.abs(geolat) <= 10.)
    assert mask.sum().item() == expected.sum().item() == 2
    assert bool((mask == expected).all())


def test_mask_dateline_crossing_box():
    grid = initialize_spherical_grid(N=6)
    # Box that genuinely straddles the ±180 antimeridian: lon 120 -> 240
    # (= -120 normalized), latitudes -20..20. Centers at xh=150 (->150) and
    # xh=210 (->-150) should be selected; rows yh in {-10, 10}.
    lons = np.array([120., 240., 240., 120., 120.])  # 240 == -120 normalized
    lats = np.array([-20., -20., 20., 20., -20.])
    mask = mask_from_grid_boundaries(lons, lats, grid)

    geolon = grid._ds["geolon"]
    geolat = grid._ds["geolat"]
    nlon = ((geolon + 180.) % 360.) - 180.
    expected = ((nlon == 150.) | (nlon == -150.)) & (np.abs(geolat) <= 10.)
    assert expected.sum().item() == 4
    assert mask.sum().item() == 4
    assert bool((mask == expected).all())
