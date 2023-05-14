import numpy as np

from .region import Region, GriddedRegion
from .boundaries import grid_boundaries_from_mask
from .overlaps import *
from .utilities import *

class Regions():
    def __init__(
        self,
        regions,
        name=None
        ):
        if type(regions) == dict:
            self.regions = regions
        else:
            raise NameError("Must provide regions to initialize.")
        if name is not None:
            self.name = name
    
    def find_all_overlaps(self, closeness_threshold=5.e3, face_indices=False):
        self.overlaps = {}
        for i, (r1name, r1) in enumerate(self.regions.items()):
            for j, (r2name, r2) in enumerate(self.regions.items()):
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
        return Regions({
            r.name: r.copy(remove_duplicate_points=remove_duplicate_points)
            for r in self.regions.values()
        })
    
class GriddedRegions(Regions):
    def __init__(
        self,
        regions,
        grid,
        name=None,
        ):
        self.grid = grid
        
        try:
            super().__init__(regions, name=name)
        except:
            raise NameError("Must provide valid regions dictionary to initialize.")
            
class MaskRegions(GriddedRegions):
    def __init__(
        self,
        mask,
        grid,
        name=None,
        ):
        self.grid = grid
        self.mask = mask
        
        i_list, j_list, lons_list, lats_list = grid_boundaries_from_mask(
            self.grid,
            mask
        )

        regions = {
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
        super().__init__(regions, grid, name=name)

def sorted_tuple(s):
    return tuple(sorted(s, key=int))