import pytest
import numpy as np
import xarray as xr
import xgcm

def initialize_spherical_grid(N=6):
    xq = np.linspace(0, 360., N+1)
    yq = np.linspace(-60, 60., N+1)
    dx = 360/N
    xh = np.linspace(0+dx/2, 360-dx/2, N)
    dy = 120/N
    yh = np.linspace(-60+dy/2, 60-dy/2, N)

    lon, lat = np.meshgrid(xh, yh)
    lon_c, lat_c = np.meshgrid(xq, yq)
    ds = xr.Dataset({}, coords={
        "xh":xr.DataArray(xh, dims=("xh",)),
        "yh":xr.DataArray(yh, dims=("yh",)),
        "xq":xr.DataArray(xq, dims=("xq",)),
        "yq":xr.DataArray(yq, dims=("yq",)),
        "geolon":xr.DataArray(lon, dims=("yh", "xh")),
        "geolat":xr.DataArray(lat, dims=("yh", "xh")),
        "geolon_c":xr.DataArray(lon_c, dims=("yq", "xq",)),
        "geolat_c":xr.DataArray(lat_c, dims=("yq", "xq",))
    })
    coords = {
        'X': {'outer': 'xq', 'center': 'xh'},
        'Y': {'outer': 'yq', 'center': 'yh'}
    }
    grid = xgcm.Grid(ds, coords=coords, padding={"X":"periodic", "Y":"extend"}, autoparse_metadata=False)
    return grid

def test_gridded_region_from_boundary():
    from regionate import GriddedRegion
    from sectionate import distance_on_unit_sphere

    lonseg = np.array([0., 120., 240, 360.])
    latseg = np.array([0.,   0.,   0.,  0.])

    grid = initialize_spherical_grid()
    region = GriddedRegion("test_region1", lonseg, latseg, grid)

    dists = distance_on_unit_sphere(
        region.lons_c,
        region.lats_c,
        np.array([0.,  60., 120., 180., 240., 300., 360.]),
        np.array([0.,   0.,   0.,   0.,   0.,   0.,   0.])
    )
    assert np.all(np.isclose(dists, 0., atol=1.e-6))

    region_rev = GriddedRegion("test_region2", lonseg[::-1], latseg[::-1], grid)
    assert np.all(np.equal(region.mask, region_rev.mask))
    
def test_curve_default_is_great_circle():
    """The default boundary->grid tracing follows great circles (geodesic edges), so
    arbitrary -- including diagonal -- boundary edges are traced correctly rather than
    collapsed onto a latitude line.

    The alternative for an axis-aligned box is `"latitude and great circle"`, which
    holds the parallel where a segment's endpoints share a latitude and takes the
    geodesic everywhere else. `"latitude circle"` on its own means every segment
    follows a parallel, so it rejects the meridional sides of a box rather than
    silently doing something else with them. On a mid-latitude box the combined
    option encloses fewer cells than the geodesic, whose edges bow poleward."""
    from regionate import GriddedRegion

    grid = initialize_spherical_grid(N=12)
    lons = np.array([60., 180., 180., 60.])   # 120-degree edges (each < 180)
    lats = np.array([-40., -40., 40., 40.])
    default = GriddedRegion("default", lons, lats, grid)
    great_circle = GriddedRegion("great", lons, lats, grid, curve="great circle")
    lat_circle = GriddedRegion("lat", lons, lats, grid,
                               curve="latitude and great circle")
    with pytest.raises(ValueError, match="constant latitude"):
        GriddedRegion("bad", lons, lats, grid, curve="latitude circle")

    assert bool((default.mask == great_circle.mask).all())  # default IS great circle
    assert bool((default.mask != lat_circle.mask).any())    # and differs from latitude circle
    # great circle bows poleward, enclosing strictly more cells than the latitude box
    assert int(great_circle.mask.sum()) > int(lat_circle.mask.sum())


def test_pole_encircling_boundary_fills_hemisphere():
    from regionate import GriddedRegion

    grid = initialize_spherical_grid()
    # A zonal line encircling the globe at the equator encloses a hemisphere; the
    # boundary cannot be drawn as a simple lon/lat polygon, so it is extended to
    # the South Pole. Both windings must enclose the same (southern) hemisphere.
    lonseg = np.array([0., 120., 240., 360.])
    latseg = np.array([0., 0., 0., 0.])
    region = GriddedRegion("hemi", lonseg, latseg, grid)
    region_rev = GriddedRegion("hemi_rev", lonseg[::-1], latseg[::-1], grid)

    assert int(region.mask.sum()) == 18              # 3 southern rows x 6 columns
    assert np.all(np.equal(region.mask, region_rev.mask))
    assert bool(region.mask.where(grid._ds.yh < 0, other=True).all())   # all southern in
    assert not bool(region.mask.where(grid._ds.yh > 0, other=False).any())  # no northern


def test_pole_cap_boundary_encloses_the_hugged_pole():
    """A pole-encircling boundary encloses the cap on the side of the pole it hugs:
    a high-northern-latitude circle encloses the NORTH cap (not its complement), a
    high-southern one the SOUTH cap. (The ambiguous equatorial case defaults to
    south -- see test_pole_encircling_boundary_fills_hemisphere.)"""
    from regionate import GriddedRegion

    grid = initialize_spherical_grid(N=6)   # row centers at -50,-30,-10,10,30,50
    lons = np.array([0., 120., 240., 360.])

    north = GriddedRegion("ncap", lons, np.full(4, 40.), grid).mask
    assert int(north.sum()) == 6                                       # only the lat=50 row
    assert bool(north.where(grid._ds.yh > 40., other=True).all())      # all northern
    assert not bool(north.where(grid._ds.yh < 40., other=False).any())  # none southern

    south = GriddedRegion("scap", lons, np.full(4, -40.), grid).mask
    assert int(south.sum()) == 6                                       # only the lat=-50 row
    assert bool(south.where(grid._ds.yh < -40., other=True).all())     # all southern
    assert not bool(south.where(grid._ds.yh > -40., other=False).any())


def test_gridded_region_from_mask():
    from regionate import MaskRegions
    
    grid = initialize_spherical_grid()
    
    # Two grid-cell wide interior region mask -> a single connected component whose
    # own boundary is one closed loop.
    mask = xr.ones_like(grid._ds.geolon).where((grid._ds.xh==90.) & (np.abs(grid._ds.yh)<=10), 0.).astype(bool)
    region_dict = MaskRegions(mask, grid).region_dict
    assert len(region_dict)==1

    region = region_dict[0]
    assert bool((region.mask == mask).all())          # component owns exactly its cells
    assert len(region.boundaries) == 1
    loop_ = region.boundaries[0]
    assert np.all(
        modequal(loop_.lons_c, np.array([ 60., 120., 120., 120.,  60.,  60.])) &
        modequal(loop_.lats_c, np.array([-20., -20.,   0.,  20.,  20.,   0.]))
    )

    # Zonal strip circling the globe: ONE connected component (it wraps the periodic-X
    # seam), whose boundary is TWO latitude circles -- the top and bottom edges, each
    # closed around the globe, not one loop cut radially at the seam. (Loop order is
    # not significant.)
    mask = xr.ones_like(grid._ds.geolon).where(np.abs(grid._ds.yh)<=10, 0.).astype(bool)
    regions = list(MaskRegions(mask, grid).region_dict.values())
    assert len(regions) == 1
    assert len(regions[0].boundaries) == 2
    assert sorted(_latitude_circle_lat(b) for b in regions[0].boundaries) == [-20., 20.]

    # Its complement is two disconnected components (a north band and a south band),
    # each an annulus with two latitude-circle loops -> four loops total at +/-20 and
    # +/-60. The +/-60 domain-wall circles survive because `_pad_center` pads a wall
    # with NaN rather than replicating the edge cell (which would make the wall look
    # like an in-mask seam and drop it).
    regions_inv = list(MaskRegions(~mask, grid).region_dict.values())
    assert len(regions_inv) == 2
    assert all(len(r.boundaries) == 2 for r in regions_inv)
    all_loops = [b for r in regions_inv for b in r.boundaries]
    assert sorted(_latitude_circle_lat(b) for b in all_loops) == [-60., -20., 20., 60.]

def test_maskregions_gr_roundtrip(tmp_path):
    """A MaskRegions collection round-trips through the multi-loop `.grs`/`.gr`
    format: each component's own mask and its full list of boundary loops (with
    their stored corner indices) survive save/reload."""
    import xgcm
    from regionate import MaskRegions
    from regionate.regions import open_grs

    grid = initialize_spherical_grid()
    # zonal strip -> ONE component with TWO boundary loops: exercises multi-loop I/O
    mask = xr.ones_like(grid._ds.geolon).where(np.abs(grid._ds.yh) <= 10, 0.).astype(bool)
    mr = MaskRegions(mask, grid, name="strip")
    mr.to_grs(f"{tmp_path}/")

    def ds_to_grid(ds):
        return xgcm.Grid(ds, coords={'X': {'outer': 'xq', 'center': 'xh'},
                                     'Y': {'outer': 'yq', 'center': 'yh'}},
                         padding={"X": "periodic", "Y": "extend"}, autoparse_metadata=False)

    reloaded = open_grs(f"{tmp_path}/strip.grs", ds_to_grid)
    assert len(reloaded.region_dict) == len(mr.region_dict) == 1
    orig = list(mr.region_dict.values())[0]
    back = list(reloaded.region_dict.values())[0]
    assert bool((back.mask.values == orig.mask.values).all())
    assert len(back.boundaries) == len(orig.boundaries) == 2

    def lat0(b):
        return round(float(np.asarray(b.lats_c)[0]), 6)
    for b in back.boundaries:
        match = next(ob for ob in orig.boundaries if lat0(ob) == lat0(b))
        assert np.array_equal(np.asarray(b.i_c), np.asarray(match.i_c))
        assert np.array_equal(np.asarray(b.j_c), np.asarray(match.j_c))
        assert np.allclose(np.mod(b.lons_c, 360.), np.mod(match.lons_c, 360.))


def modequal(a,b):
    return np.equal(np.mod(a, 360.), np.mod(b, 360.))

def _latitude_circle_lat(region):
    """Assert `region`'s boundary loop circles the globe at a single latitude, spanning
    every longitude once (periodic-X seam stitched, no radial cut), and return that lat."""
    lats = np.asarray(region.lats_c)
    lons = set(np.mod(np.round(np.asarray(region.lons_c)), 360.).tolist())
    assert np.allclose(lats, lats[0]), lats
    assert lons == {0., 60., 120., 180., 240., 300.}, lons
    return float(lats[0])