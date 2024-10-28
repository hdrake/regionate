import geopandas as gpd
from shapely.geometry import Polygon
import regionmask

import sectionate as sec
import numpy as np

from .utilities import *

def get_region_boundary_grid_indices(lons, lats, grid):
    """Find boundary coordinates and grid indices that approximate a polygon.

    ARGUMENTS
    ---------
    lons : list or np.ndarray of longitudes
    lats : list or np.ndarray of latitudes
    grid : `xgcm.Grid` instance
    
    RETURNS
    -------
    (i, j, lons_c, lats_c, lons_uv, lats_uv)

    i : "X"-axis grid indices of corner points
    j : "Y"-axis grid indices of corner points
    lons_c : longitudes of corner points
    lats_c : latitudes of corner points
    lons_uv : longitudes of (u,v) velocity faces
    lats_uv : latitudes of (u,v) velocity faces
    """
    if (lons[0], lats[0]) != (lons[-1], lats[-1]):
        lons, lats = loop(lons), loop(lats)
        
    i, j, lons_c, lats_c = sec.grid_section(
        grid,
        lons,
        lats,
        topology="MOM-tripolar"
    )
    
    uvindices = sec.uvindices_from_qindices(grid, i, j)
    
    lons_uv, lats_uv = sec.uvcoords_from_uvindices(
        grid,
        uvindices
    )

    return (i, j, lons_c[:-1], lats_c[:-1], lons_uv, lats_uv)

def mask_from_grid_boundaries(
    lons_c,
    lats_c,
    grid,
    along_boundary=False,
    coordnames={'h': ('geolon', 'geolat')},
    ):
    """Find mask bounded by a sequence of cell corner coordinates

    ARGUMENTS
    ---------
    lons_c [list or np.ndarray] -- cell corner longitudes
    lats_c [list or np.ndarray] -- cell corner latitudes
    grid [xgcm.Grid] -- ocean model grid
    
    RETURNS
    -------
    region_grid_mask : np.ndarray of bool type
    """
    
    Δlon = np.sum(np.diff(lons_c)[np.abs(np.diff(lons_c)) < 180])
    
    if along_boundary:
        polygon_geom = Polygon(zip(lons_c, lats_c))

        crs = 'epsg:4326'
        polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon_geom])
        region_grid_mask = ~np.isnan(
            regionmask.mask_geopandas(
                polygon,
                grid._ds[coordnames['h'][0]],
                lat=grid._ds[coordnames['h'][1]],
                wrap_lon=False
            )
        )
        
    elif np.abs(Δlon) < 180.:
        lons = loop(lons_c)
        lats = loop(lats_c)
        wrapped_lons = wrap_continuously(lons)
        minlon = np.min(wrapped_lons)
        polygon_geom = Polygon(zip(np.mod(wrapped_lons-minlon, 360.), lats))

        crs = 'epsg:4326'
        polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon_geom])
        region_grid_mask = ~np.isnan(
            regionmask.mask_geopandas(
                polygon,
                np.mod(grid._ds[coordnames['h'][0]]-minlon, 360.),
                lat=grid._ds[coordnames['h'][1]],
                wrap_lon='360'
            )
        )
        
    else:
        # Follow sectionate convention for orientations of polygons
        # on stereographic plane (relative to South Pole)
        s = np.sign(Δlon).astype(int)
        if s==-1:
            lons_c = lons_c[::-1]
            lats_c = lats_c[::-1]
            s = 1

        min_idx = np.argmin(lons_c)

        lons = np.roll(lons_c, -min_idx)
        lats = np.roll(lats_c, -min_idx)

        lons = np.append(lon_mod(lons[-1], lons[0]), lons)
        lats = np.append(lats[-1], lats)

        diffs = s*(lons[np.newaxis, :] - lons[:, np.newaxis])
        diffs[np.tril_indices(lons.size)]*=-1
        single_valued = ~np.any(diffs < 0, axis=1)

        roll_idx = np.argmax(single_valued[::s])
        lons = np.roll(lons[::s], -roll_idx)[::s]
        lats = np.roll(lats[::s], -roll_idx)[::s]
        lons[::s][-roll_idx:] = lons[::s][-roll_idx:]

        # Make sure that we get everything to the South
        min_idx = np.argmin(lons)
        max_idx = np.argmax(lons)
        lons = np.append(
            lons, [
                lons[max_idx]+10,
                lons[max_idx]+10,
                lons[min_idx]-10,
                lons[min_idx]-10
            ]
        )
        lats = np.append(
            lats, [
                lats[max_idx],
                -90,
                -90,
                lats[min_idx]
            ]
        )
        
        polygon_geom = Polygon(zip(lons, lats))
        crs = 'epsg:4326'
        polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon_geom])
        region_grid_mask = ~np.isnan(
            regionmask.mask_geopandas(
                polygon,
                grid._ds[coordnames['h'][0]],
                lat=grid._ds[coordnames['h'][1]],
                wrap_lon=False
            )
        )
    
    return region_grid_mask

def wrap_continuously(x, limit_discontinuity=180.):
    new_x = x.copy()
    for i in range(len(new_x)-1):
        if new_x[i+1]-new_x[i] >= 180.:
            new_x[i+1] -= 360.
        elif new_x[i+1]-new_x[i] < -180:
            new_x[i+1] += 360.
    return new_x
