import contourpy
import xarray as xr
import numpy as np

from .utilities import loop, unique_list
from sectionate.transports import check_symmetric

def grid_boundaries_from_mask(grid, mask):
    """Use `contourpy` to get indices of cell corners that bound the mask

    Returns lists with a common length that is the number of discrete contours
    that bound the mask.
    
    ARGUMENTS
    ---------
    grid : `xgcm.Grid` instance
    mask : `xr.DataArray` instance of type `bool`
    
    RETURNS
    -------
    i_c_list, j_c_list, lons_c_list, lats_c_list
    """
    
    symmetric = check_symmetric(grid)
    
    contours = (
        contourpy.contour_generator(
            np.arange(-1, mask.shape[1]+1),
            np.arange(-1, mask.shape[0]+1),
            np.pad(mask.values, np.array([1,1]))
        )
        .create_contour(0.5)
    )
    i_c_list = []
    j_c_list = []
    lons_c_list = []
    lats_c_list = []
    for c in contours:
        i_c, j_c = c[:-1,0], c[:-1,1]
        
        i_c_new, j_c_new = i_c.copy(), j_c.copy()
        
        i_inc = np.roll(i_c, -1)-i_c
        j_inc = np.roll(j_c, -1)-j_c
        
        i_c_new[(i_c%1)==0.0] = (i_c - (i_inc<0))[(i_c%1)==0.0] + symmetric
        j_c_new[(j_c%1)==0.0] = (j_c - (j_inc<0))[(j_c%1)==0.0] + symmetric
        i_c_new[(i_c%1)==0.5] = np.floor(i_c[(i_c%1)==0.5]) + symmetric
        j_c_new[(j_c%1)==0.5] = np.floor(j_c[(j_c%1)==0.5]) + symmetric
        
        i_c_new, j_c_new = loop(i_c_new).astype(np.int64), loop(j_c_new).astype(np.int64)
                        
        i_c_list.append(i_c_new)
        j_c_list.append(j_c_new)
        
        idx = {
            "xq":xr.DataArray(i_c_new, dims=("pt",)),
            "yq":xr.DataArray(j_c_new, dims=("pt",))
        }
        lons_c_list.append(grid._ds["geolon_c"].isel(idx).values[:-1])
        lats_c_list.append(grid._ds["geolat_c"].isel(idx).values[:-1])
        
    return i_c_list, j_c_list, lons_c_list, lats_c_list