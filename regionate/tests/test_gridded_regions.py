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
    grid = xgcm.Grid(ds, coords=coords, boundary={"X":"periodic", "Y":"extend"}, autoparse_metadata=False)
    return grid

def test_gridded_region_from_boundary():
    from regionate import GriddedRegion
    from sectionate import distance_on_unit_sphere

    lonseg = np.array([0., 120., 240, 360.])
    latseg = np.array([0., 0., 0., 0.])

    grid = initialize_spherical_grid()
    region = GriddedRegion("test_region1", lonseg, latseg, grid)

    dists = distance_on_unit_sphere(
        region.lons_c,
        region.lats_c,
        np.array([360.,  60., 120., 180., 240., 300.]),
        np.array([0., 0., 0., 0., 0., 0.])
    )
    assert np.all(np.isclose(dists, 0., atol=1.e-6))

    region_rev = GriddedRegion("test_region2", lonseg[::-1], latseg[::-1], grid)
    assert np.all(np.equal(region.mask, region_rev.mask))
    
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
    
    # Zonal strip mask
    mask = xr.ones_like(grid._ds.geolon).where(np.abs(grid._ds.yh)<=10, 0.).astype(bool)
    region = MaskRegions(mask, grid).region_dict[0]
    assert np.all(
        modequal(region.lons_c, np.array([360.,  60., 120., 180., 240., 300., 360., 360., 360., 300., 240., 180., 120.,  60., 360., 360.])) &
        modequal(region.lats_c, np.array([-20., -20., -20., -20., -20., -20., -20.,   0.,  20.,  20.,  20., 20.,  20.,  20.,  20.,   0.]))
    )
    
    # All but zonal strip mask (two separate regions outside)
    region_inv = MaskRegions(~mask, grid).region_dict
    assert np.all(
        modequal(region_inv[0].lons_c, np.array([360.,  60., 120., 180., 240., 300., 360., 360., 360., 300., 240., 180., 120.,  60., 360., 360.])) &
        modequal(region_inv[0].lats_c, np.array([-60., -60., -60., -60., -60., -60., -60., -40., -20., -20., -20., -20., -20., -20., -20., -40.])) &
        modequal(region_inv[1].lons_c, np.array([360.,  60., 120., 180., 240., 300., 360., 360., 360., 300., 240., 180., 120.,  60., 360., 360.])) &
        modequal(region_inv[1].lats_c, np.array([20.,   20.,  20.,  20.,  20.,  20.,  20.,  40.,  60.,  60.,  60.,  60.,  60.,  60.,  60.,  40.]))
    )

def modequal(a,b):
    return np.equal(np.mod(a, 360.), np.mod(b, 360.))