import numpy as np
import xarray as xr

from .region import Region, GriddedRegion, BoundedRegion, open_gr
from .boundaries import grid_boundaries_from_mask
from .overlaps import *
from .utilities import *

import os
from pathlib import Path

class Regions():
    """
    A dictionary of polygonal regions defined by a list or array of geographical coordinates.
    """
    def __init__(
        self,
        region_dict,
        name=None
        ):
        """
        Create a `Regions` object from a dictionary mapping region names to `Region` instances.

        PARAMETERS
        ----------
        region_dict : dictionary mapping region `name` (str) to `Region` instance
        name : str or None (default: None)
            Overarching name of the collection of regions

        RETURNS
        -------
        `Regions` instance

        Examples
        --------
        >>> lons, lats = np.array([-80., -66., -65.]), np.array([ 26.,  18.,  32.])
        >>> region = reg.Region('Bermuda Triangle', lons, lats)
        >>> regions = reg.Regions({region.name: region})
        """
        if type(region_dict) == dict:
            self.region_dict = region_dict
        else:
            raise NameError("Must provide `regions_dict` to initialize.")
        if name is not None:
            self.name = name
    
    def find_all_overlaps(self, closeness_threshold=5.e3, face_indices=False):
        """Finds indices of `self.overlaps` between region boundaries

        Combined with `sectionate`, this can be used to determine transports across
        shared boundaries between two adjacent regions.
        """
        self.overlaps = {}
        for i, (r1name, r1) in enumerate(self.region_dict.items()):
            for j, (r2name, r2) in enumerate(self.region_dict.items()):
                if r1name<r2name:
                    overlaps = find_indices_of_overlaps(
                        r1,
                        r2,
                        closeness_threshold=closeness_threshold,
                        face_indices=face_indices
                    )
                    if len(overlaps[r1name]):
                        self.overlaps[sorted_tuple((r1name, r2name))] = overlaps
                
    def copy(self, remove_duplicate_points=False):
        """
        Returns a copy of the Regions.

        PARAMETERS
        ----------
        remove_duplicate_points : bool
            Default: False. If True, prunes any duplicate points from the input arrays (lons, lats).
            
        RETURNS
        ----------
        region_copy : `regionate.regions.Regions` type
            Copy of the region.
        """
        return Regions({
            r.name: r.copy(remove_duplicate_points=remove_duplicate_points)
            for r in self.region_dict.values()
        })
    
class GriddedRegions(Regions):
    """
    A dictionary of named polygonal regions that exactly conform to the velocity faces of a C-grid ocean model.
    """
    def __init__(
        self,
        region_dict,
        grid,
        name=None,
        ):
        """
        Create a `GriddedRegions` object from a dictionary mapping region names to `GriddedRegion` instances.

        PARAMETERS
        ----------
        region_dict : dictionary mapping region `name` (str) to `Region`, `GriddedRegion`, or `BoundedRegion` instance
        name : str or None (default: None)
            Overarching name of the collection of regions

        RETURNS
        -------
        `GriddedRegions` instance
        """
        self.grid = grid
        
        super().__init__(region_dict, name=name)
        try:
            if all([type(v) in [Region, GriddedRegion, BoundedRegion] for v in region_dict.values()]):
                super().__init__(region_dict, name=name)
            else:
                raise NameError("""Values in `region_dict` dictionary must be instances of
                `Region`, `GriddedRegion`, or `BoundedRegion`.""")
        except:
            raise NameError("Must provide valid `region_dict` dictionary to initialize.")

    def to_grs(self, path):
        """
        TO DO
        """
        
        # Create .grs file directory
        grs_path = f"{path}{self.name.replace(' ','_')}.grs"
        Path(grs_path).mkdir(parents=True, exist_ok=True)

        # Write grid dataset (without variables) to NetCDF file
        grid = self.grid
        grid_path = f"{grs_path}/grid.nc"
        grid._ds.drop_vars([v for v in grid._ds.data_vars]).to_netcdf(grid_path)

        # Write boundary information for each region to separate .gr file directory
        for region_name, region in self.region_dict.items():
            region.to_gr(f"{grs_path}/")

class MaskRegions(GriddedRegions):
    """
    A dictionary of polygonal regions that exactly conform to the velocity faces bounding a mask in a C-grid ocean model.
    """
    def __init__(
        self,
        mask,
        grid,
        name=None,
        ):
        """
        Create a `MaskRegions` object from a mask and accompanying `xgcm.Grid` instance.

        PARAMETERS
        ----------
        mask : None or xr.DataArray (default: None)
            If None, does not apply any mask.
        grid : `xgcm.Grid` instance
        name : str or None (default: None)
            Overarching name of the collection of regions

        RETURNS
        -------
        `GriddedRegions` instance
        """
        
        self.grid = grid
        self.mask = mask
        
        i_list, j_list, lons_list, lats_list = grid_boundaries_from_mask(
            self.grid,
            mask
        )

        region_dict = {
            r_num: GriddedRegion(
                str(r_num),
                lons,
                lats,
                self.grid,
                mask=mask,
                ij=(i,j)
            )
            for r_num, (i, j, lons, lats)
            in enumerate(zip(i_list, j_list, lons_list, lats_list))
        }
        super().__init__(region_dict, grid, name=name)

def open_grs(path, ds_to_grid):
    """
    Creates a GriddedRegions instance from a .grs file format

    PARAMETERS
    ----------
    path : str
        path to .grs directory
    ds_to_grid : function
        a function whose only argument is an xr.Dataset object.
        This function turns `xr.open_dataset(f"{path}/grid.nc")`
        into a corresponding `xgcm.Grid` object.

    RETURNS
    -------
    `GriddedRegions` instance
    """
    name = path.split('/')[-1][:-4]
    grid = ds_to_grid(xr.open_dataset(f"{path}/grid.nc"))
    grs_files = [f for f in os.listdir(f"{path}/") if f!='grid.nc']
    region_dict = {}
    for gr_file in grs_files:
        gr_name = gr_file[:-3]
        region_dict[gr_name] = open_gr(f"{path}/{gr_file}", ds_to_grid)

    return GriddedRegions(region_dict, grid, name=name)

def sorted_tuple(s):
    """Sort tuples (by integer)"""
    return tuple(sorted(s, key=int))