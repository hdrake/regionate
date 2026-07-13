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


def test_gridded_region_from_mask():
    from regionate import MaskRegions
    
    grid = initialize_spherical_grid()
    
    # Two grid-cell wide interior region mask
    mask = xr.ones_like(grid._ds.geolon).where((grid._ds.xh==90.) & (np.abs(grid._ds.yh)<=10), 0.).astype(bool)
    region_dict = MaskRegions(mask, grid).region_dict
    assert len(region_dict)==1
    
    region = region_dict[0]
    assert np.all(
        modequal(region.lons_c, np.array([ 60., 120., 120., 120.,  60.,  60.])) &
        modequal(region.lats_c, np.array([-20., -20.,   0.,  20.,  20.,   0.]))
    )
    
    # Zonal strip circling the globe: because the tracer stitches across the periodic-X
    # seam, its boundary is TWO latitude circles (the top and bottom edges), each closed
    # around the globe -- not one loop cut radially at the seam. (Loop order is not
    # significant.)
    mask = xr.ones_like(grid._ds.geolon).where(np.abs(grid._ds.yh)<=10, 0.).astype(bool)
    regions = list(MaskRegions(mask, grid).region_dict.values())
    assert len(regions) == 2
    assert sorted(_latitude_circle_lat(r) for r in regions) == [-20., 20.]

    # Its complement is two disconnected bands (north and south), each an annulus -> four
    # latitude circles at +/-20 and +/-60. The +/-60 domain-wall circles survive because
    # `_pad_center` pads a wall with NaN rather than replicating the edge cell (which
    # would make the wall look like an in-mask seam and drop it).
    regions_inv = list(MaskRegions(~mask, grid).region_dict.values())
    assert len(regions_inv) == 4
    assert sorted(_latitude_circle_lat(r) for r in regions_inv) == [-60., -20., 20., 60.]

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