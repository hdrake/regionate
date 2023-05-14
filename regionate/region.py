import numpy as np

from sectionate import is_section_counterclockwise
from .utilities import *
from .grid_conform import (
    get_region_boundary_grid_indices,
    mask_from_grid_boundaries
)

class Region:
    def __init__(
        self,
        name,
        lons,
        lats,
        counterclockwise=None,
        force_ccw=False,
        remove_duplicate_points=False
        ):
        self.name = name
        self.lons = lons
        self.lats = lats
        
        if remove_duplicate_points:
            self.remove_duplicate_points()

        if counterclockwise is not None:
            self.counterclockwise = counterclockwise
        else:
            self.counterclockwise = is_section_counterclockwise(
                loop(self.lons),
                loop(self.lats)
            )
            
        if force_ccw:
            self.make_counterclockwise()

    def copy(self, remove_duplicate_points=False):
        return Region(
            self.name,
            self.lons.copy(),
            self.lats.copy(),
            counterclockwise=self.counterclockwise,
            remove_duplicate_points=remove_duplicate_points
        )
    
    def make_counterclockwise(self):
        if not(self.counterclockwise):
            self.lons = self.lons[::-1]
            self.lats = self.lats[::-1]
            self.counterclockwise = True

    def remove_duplicate_points(self, closeness_threshold=5.e3):
        self.lons, self.lats = unique_lonlat(self.lons, self.lats, closeness_threshold=closeness_threshold)
        
class GriddedRegion(Region):
    def __init__(
        self,
        name,
        lons,
        lats,
        grid,
        positive_in=True,
        mask=None,
        ij=None
        ):
        self.grid = grid
        
        if len(lons)>=3 and len(lats)>=3 and ij is None:
            self.initiate_from_boundary(lons, lats, mask=mask, positive_in=positive_in)
        elif ij is None:
            raise NameError("Must provide lons and lats as lists or arrays\
            to define the region.")
        else:
            self.lons = lons
            self.lats = lats
            self.i = ij[0]
            self.j = ij[1]
            self.mask = mask
        
        super().__init__(
            name=name,
            lons=self.lons,
            lats=self.lats
        )
        
    def initiate_from_boundary(
        self,
        lons,
        lats,
        positive_in=True,
        mask=None,
        ):

        self.i, self.j, self.lons, self.lats, self.lons_uv, self.lats_uv = (
            get_region_boundary_grid_indices(
                lons.copy(),
                lats.copy(),
                self.grid,
            )
        )
        if mask is None:
            mask = mask_from_grid_boundaries(
                self.lons,
                self.lats,
                self.grid
            )
        self.mask = mask.astype(bool) ^ (not positive_in)