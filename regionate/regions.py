import numpy as np

from .region import Region, GriddedRegion
from .boundaries import grid_boundaries_from_mask
from .overlaps import *
from .utilities import *

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
        >>> region = reg.Region("Bermuda Triangle", lons, lats)
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
        region_dict : dictionary mapping region `name` (str) to `Region` or `GriddedRegion` instance
        name : str or None (default: None)
            Overarching name of the collection of regions

        RETURNS
        -------
        `GriddedRegions` instance
        """
        self.grid = grid
        
        try:
            if all([type(v) in [Region, GriddedRegion] for v in region_dict.values()]):
                super().__init__(region_dict, name=name)
            else:
                raise NameError("Values in `region_dict` dictionary must be of instances of `Region` or `GriddedRegion`.")
        except:
            raise NameError("Must provide valid `region_dict` dictionary to initialize.")
            
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

def sorted_tuple(s):
    """Sort tuples (by integer)"""
    return tuple(sorted(s, key=int))