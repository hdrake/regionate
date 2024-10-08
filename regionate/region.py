import numpy as np

from sectionate import is_section_counterclockwise
from .utilities import *
from .grid_conform import (
    get_region_boundary_grid_indices,
    mask_from_grid_boundaries
)

class Region:
    """
    A named polygonal region defined by a list or array of geographical coordinates.
    """
    def __init__(
        self,
        name,
        lons,
        lats,
        force_ccw=False,
        remove_duplicate_points=False
        ):
        """
        Create a Region object (named `name`) from arrays of (`lons`, `lats`).

        PARAMETERS
        ----------
        name : str
            Name of the region
        lons : list or np.ndarray
            Longitudes (in degrees).
        lats : list or np.ndarray
            Latitudes (in degrees).
        force_ccw : bool
            Default: False. If True, checks if Region is clockwise and, if it is,
            swaps the order of the points that define it such that it becomes counterclockwise.
        remove_duplicate_points : bool
            Default: False. If True, prunes any duplicate points from the input arrays (lons, lats).

        RETURNS
        -------
        Region instance

        Examples
        --------
        >>> lons, lats = np.array([-80., -66., -65.]), np.array([ 26.,  18.,  32.])
        >>> region = reg.Region("Bermuda Triangle", lons, lats)
        """
        
        self.name = name
        self.lons = lons
        self.lats = lats
        
        if remove_duplicate_points:
            self.remove_duplicate_points()

        self.counterclockwise = is_section_counterclockwise(
            loop(self.lons),
            loop(self.lats),
            geometry="spherical"
        )
            
        if force_ccw:
            self.make_counterclockwise()

    def copy(self, remove_duplicate_points=False):
        """
        Returns a copy of the Region.

        PARAMETERS
        ----------
        remove_duplicate_points : bool
            Default: False. If True, prunes any duplicate points from the input arrays (lons, lats).
            
        RETURNS
        ----------
        region_copy : regionate.region.Region
            Copy of the region.
        """
        return Region(
            self.name,
            self.lons.copy(),
            self.lats.copy(),
            remove_duplicate_points=remove_duplicate_points
        )
    
    def make_counterclockwise(self):
        """
        Checks if the section is clockwise and flips its direction if it is to make it counterclockwise.
        """
        if not(self.counterclockwise):
            self.lons = self.lons[::-1]
            self.lats = self.lats[::-1]
            self.counterclockwise = True

    def remove_duplicate_points(self, closeness_threshold=5.e3):
        """
        Removes any duplicate points.
        
        PARAMETERS
        ----------
        closeness_threshold : float
            A short distance within which points are deemed to be identical. Default: 5.e3.
        """
        self.lons, self.lats = unique_lonlat(self.lons, self.lats, closeness_threshold=closeness_threshold)
        
class GriddedRegion(Region):
    """
    A named polygonal region that exactly conforms to the velocity faces of a C-grid ocean model.
    """
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
        """
        Create a Region object (named `name`) from arrays of (`lons`, `lats`) and an ocean model `grid`. 

        PARAMETERS
        ----------
        name : str
            Name of the region
        lons : list or np.ndarray
            Longitudes (in degrees).
        lats : list or np.ndarray
            Latitudes (in degrees).
        grid : xgcm.Grid
        positive_in : bool
            Default: True. If True, prunes any duplicate points from the input arrays (lons, lats).
        mask : None or xr.DataArray (default: None)
            If None, does not apply any mask.
        ij : None or list
            If None, the indices of grid coordinates closest to provided coordinates `self.i` and `self.j`
            are inferred from the model grid. If list, assume two elements in the list and
            extract `self.i = ij[0]` and `self.j = ij[1]`.

        RETURNS
        -------
        GriddedRegion instance

        Examples
        --------
        >>> grid = xgcm.Grid(...) # TO DO: minimal example grid
        >>> lons, lats = np.array([-80., -66., -65.]), np.array([ 26.,  18.,  32.])
        >>> region = reg.Region("Bermuda Triangle", lons, lats, grid)
        """
        self.grid = grid
        
        if len(lons)>=3 and len(lats)>=3 and ij is None:
            self.initiate_from_boundary(
                lons,
                lats,
                mask=mask,
                positive_in=positive_in
            )
        elif ij is None:
            raise NameError("Must provide lons and lats as lists or arrays\
            to define the region.")
        else:
            self.lons = lons
            self.lats = lats
            self.i = ij[0]
            self.j = ij[1]
            if mask is None:
                self.mask = mask_from_grid_boundaries(
                    self.lons,
                    self.lats,
                    self.grid,
                    along_boundary=True
                )
            else:
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
        """
        TO DO
        """

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