import contourpy
import xarray as xr
import numpy as np

from .utilities import loop, unique_list
from sectionate.transports import check_symmetric

def grid_boundaries_from_mask(grid, mask):
    
    symmetric = check_symmetric(grid)
    
    contours = (
        contourpy.contour_generator(
            np.arange(-1, mask.shape[1]+1),
            np.arange(-1, mask.shape[0]+1),
            np.pad(mask.values, np.array([1,1]))
        )
        .create_contour(0.5)
    )
    i_list = []
    j_list = []
    lons_list = []
    lats_list = []
    for c in contours:
        i, j = c[:-1,0], c[:-1,1]
        
        i_new, j_new = i.copy(), j.copy()
        
        i_inc = np.roll(i, -1)-i
        j_inc = np.roll(j, -1)-j
        
        i_new[(i%1)==0.0] = (i - (i_inc<0))[(i%1)==0.0] + symmetric
        j_new[(j%1)==0.0] = (j - (j_inc<0))[(j%1)==0.0] + symmetric
        i_new[(i%1)==0.5] = np.floor(i[(i%1)==0.5]) + symmetric
        j_new[(j%1)==0.5] = np.floor(j[(j%1)==0.5]) + symmetric
        
        i_new, j_new = loop(i_new).astype(np.int64), loop(j_new).astype(np.int64)
                        
        i_list.append(i_new)
        j_list.append(j_new)
        
        idx = {
            "xq":xr.DataArray(i_new, dims=("pt",)),
            "yq":xr.DataArray(j_new, dims=("pt",))
        }
        lons_list.append(grid._ds["geolon_c"].isel(idx).values[:-1])
        lats_list.append(grid._ds["geolat_c"].isel(idx).values[:-1])
        
    return i_list, j_list, lons_list, lats_list