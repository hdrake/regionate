import geopandas as gpd
from shapely.geometry import Polygon
import regionmask

import sectionate as sec
import numpy as np

def conform_basin_to_ocean_grid(b, ocean_grid):
    b.lons_input, b.lats_input = b.lons.copy(), b.lats.copy()
    b.i, b.j, b.lons, b.lats, b.lons_uv, b.lats_uv = basin_boundary_grid_indices(b, ocean_grid)
    b.mask = basin_interior_grid_mask(b, ocean_grid)
    return b

def basin_boundary_grid_indices(
    b, ocean_grid,
    coordnames={
        'h': ('geolon',   'geolat'  ),
        'q': ('geolon_c', 'geolat_c'),
    }):
    
    symmetric = ocean_grid[coordnames['h'][0]].shape!=ocean_grid[coordnames['q'][0]].shape
    
    i, j, lons_c, lats_c = sec.create_section_composite(
        ocean_grid[coordnames['q'][0]],
        ocean_grid[coordnames['q'][1]],
        np.append(b.lons, b.lons[0]),
        np.append(b.lats, b.lats[0]),
        symmetric,
        closed=True
    )
    
    uvindices = sec.uvindices_from_qindices(i, j, symmetric)
    
    lons_uv, lats_uv = sec.uvcoords_from_uvindices(
        ocean_grid,
        uvindices,
    )

    return i, j, lons_c[:-1], lats_c[:-1], lons_uv, lats_uv

def basin_interior_grid_mask(
    b, ocean_grid,
    coordnames={
        'h': ('geolon', 'geolat'),
    }):
    Δlon = np.sum(np.diff(b.lons)[np.abs(np.diff(b.lons)) < 180])
    
    if np.abs(Δlon) < 180.:
        lons = np.append(b.lons, [b.lons[0]])
        lats = np.append(b.lats, [b.lats[0]])
        wrapped_lons = wrap_continuously(lons)
        minlon = np.min(wrapped_lons)
        polygon_geom = Polygon(zip(np.mod(wrapped_lons-minlon, 360.), lats))

        crs = 'epsg:4326'
        polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon_geom])
        basin_grid_mask = ~np.isnan(regionmask.mask_geopandas(polygon, np.mod(ocean_grid[coordnames['h'][0]]-minlon, 360.), lat=ocean_grid[coordnames['h'][1]], wrap_lon='360'))
        
    else:
        # Flag specifies whether circumgeo regions are bounded by north or south pole
        s = np.sign(Δlon).astype(int)
        
        if s>0:
            min_idx = np.argmin(b.lons)
        else:
            min_idx = np.argmax(b.lons)
        lons = np.roll(b.lons, -min_idx)
        lats = np.roll(b.lats, -min_idx)

        diffs = s*(lons[np.newaxis, :] - lons[:, np.newaxis])
        diffs[np.tril_indices(lons.size)]*=-1
        single_valued = ~np.any(diffs < 0, axis=1)
        
        roll_idx = np.argmax(single_valued[::s])
        lons = np.roll(lons[::s], -roll_idx)[::s]
        lats = np.roll(lats[::s], -roll_idx)[::s]
        lons[::s][-roll_idx:] = lons[::s][-roll_idx:] - s*360.

        min_idx = np.argmin(lons)
        lons = np.append(lons, [lons[-1], lons[-1]-s*360, lons[-1]-s*360, lons[0]])
        lats = np.append(lats, [s*90,     s*90,           lats[0],        lats[0]])

        minlon = np.min(lons)
        polygon_geom = Polygon(zip(lons-minlon, lats))
        crs = 'epsg:4326'
        polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon_geom])
        basin_grid_mask = ~np.isnan(regionmask.mask_geopandas(polygon, ocean_grid[coordnames['h'][0]]-minlon, lat=ocean_grid[coordnames['h'][1]], wrap_lon='360'))
    
    return basin_grid_mask

def wrap_continuously(x, limit_discontinuity=180.):
    new_x = x.copy()
    for i in range(len(new_x)-1):
        if new_x[i+1]-new_x[i] >= 180.:
            new_x[i+1] -= 360.
        elif new_x[i+1]-new_x[i] < -180:
            new_x[i+1] += 360.
    return new_x