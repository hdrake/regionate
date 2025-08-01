import numpy as np
import xarray as xr

import sectionate as sec
from sectionate import is_section_counterclockwise
from .utilities import *
from .grid_conform import (
    get_region_boundary_grid_indices,
    mask_from_grid_boundaries
)

import os
from pathlib import Path

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
        >>> region = reg.Region('Bermuda Triangle', lons, lats)
        """
        
        self.name = name
        self.lons_c = lons
        self.lats_c = lats
        
        if remove_duplicate_points:
            self.remove_duplicate_points()

        self.counterclockwise = is_section_counterclockwise(
            loop(self.lons_c),
            loop(self.lats_c),
            geometry='spherical'
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
            self.lons_c.copy(),
            self.lats_c.copy(),
            remove_duplicate_points=remove_duplicate_points
        )
    
    def make_counterclockwise(self):
        """
        Checks if the section is clockwise and flips its direction if it is to make it counterclockwise.
        """
        if not(self.counterclockwise):
            self.lons_c = self.lons_c[::-1]
            self.lats_c = self.lats_c[::-1]
            self.counterclockwise = True

    def remove_duplicate_points(self, closeness_threshold=5.e3):
        """
        Removes any duplicate points.
        
        PARAMETERS
        ----------
        closeness_threshold : float
            A short distance within which points are deemed to be identical. Default: 5.e3.
        """
        self.lons_c, self.lats_c = unique_lonlat(
            self.lons_c,
            self.lats_c,
            closeness_threshold=closeness_threshold
        )

    def __repr__(self):
        return f"{str(type(self))[8:-2]}('{self.name}')"
        
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
        ij=None,
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
            If None, the indices of grid coordinates closest to provided coordinates `self.i_c` and `self.j_c`
            are inferred from the model grid. If list, assume two elements in the list and
            extract `self.i_c = ij[0]` and `self.j_c = ij[1]`.

        RETURNS
        -------
        GriddedRegion instance

        Examples
        --------
        >>> grid = xgcm.Grid(...) # TO DO: minimal example grid
        >>> lons, lats = np.array([-80., -66., -65.]), np.array([ 26.,  18.,  32.])
        >>> region = reg.Region('Bermuda Triangle', lons, lats, grid)
        """
        self.grid = grid
        self.save = {}
        
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
            self.lons_c = lons
            self.lats_c = lats
            self.i_c = ij[0]
            self.j_c = ij[1]
            if mask is None:
                self.mask = mask_from_grid_boundaries(
                    self.lons_c,
                    self.lats_c,
                    self.grid,
                    along_boundary=True
                )
            else:
                self.mask = mask
        
        super().__init__(
            name=name,
            lons=self.lons_c,
            lats=self.lats_c
        )
        
    def initiate_from_boundary(
        self,
        lons,
        lats,
        positive_in=True,
        mask=None
        ):
        """
        TO DO
        """

        self.i_c, self.j_c, self.lons_c, self.lats_c, self.lons_uv, self.lats_uv = (
            get_region_boundary_grid_indices(
                lons.copy(),
                lats.copy(),
                self.grid
            )
        )
        if mask is None:
            mask = mask_from_grid_boundaries(
                self.lons_c,
                self.lats_c,
                self.grid
            )
        self.mask = mask.astype(bool) ^ (not positive_in)

    def to_gr(self, path):
        """Save the GriddedRegion object in a .gr format file directory

        There are two key files within each .gr directory:
          - a `grid.nc` file that contains information about the coordinates
          requires to create an `xgcm.Grid` instance
          - a `boundary.nc` file that contains the region's tracer cell mask
          and the coordinates and indices of the corner cells that define its
          boundary.

        To do:
          - Subdirectory for child boundary information

        Arguments
        ---------
        path [str] -- path to directory where the .gr file directory
            should be saved. The filename will be [GriddedRegion.name].gr

        Example
        -------
        >>> gridded_region.save('../data/')
        """
        gr_path = f"{path}/{self.name.replace(' ','_')}.gr/"
        Path(gr_path).mkdir(parents=True, exist_ok=True)

        grid_path = f"{gr_path}/grid.nc"
        parent_grid = f"{path}/../grid.nc"
        if os.path.isfile(parent_grid):
            os.symlink(parent_grid, grid_path)
        else:
            grid = self.grid
            grid._ds.drop_vars([v for v in grid._ds.data_vars]).to_netcdf(grid_path)

        # Write boundary information to NetCDF file
        vertex = xr.DataArray(np.arange(0, self.i_c.size), dims=('vertex',))
        face = xr.DataArray(np.arange(0.5, self.i_c.size-1), dims=('face',))
        ds = xr.Dataset({}, coords={'vertex': vertex, 'face': face})
        for v in ['lons', 'lats', 'i', 'j']:
            var = getattr(self,v)
            if v in ['lons', 'lats']:
                var = loop(var)
            ds[v] = xr.DataArray(var, dims=('vertex',))
        for v in ['lons_uv', 'lats_uv']:
            ds[v] = xr.DataArray(getattr(self,v), dims=('face',))
        ds['mask'] = self.mask
        ds.to_netcdf(f"{gr_path}/region.nc")

        for (k,v) in self.save.items():
            v.to_netcdf(f"{gr_path}/{k}.nc")

        # Write boundary information for each child section to NetCDF file
        child_path = f"{gr_path}/children/"
        Path(child_path).mkdir(parents=True, exist_ok=True)
        for child in self.children.values():
            child_name = child.name.replace(' ','_')
            sec_path = f"{gr_path}/children/{child_name}.sec"
            Path(sec_path).mkdir(parents=True, exist_ok=True)
            
            vertex = xr.DataArray(np.arange(0, child.i_c.size), dims=('vertex',))
            ds = xr.Dataset({}, coords={'vertex': vertex})
            for k in ['lons', 'lats', 'i', 'j']:
                var = getattr(child,k)
                ds[k] = xr.DataArray(var, dims=('vertex',))
            ds.to_netcdf(f"{sec_path}/section.nc")
            for k, ds_save in child.save.items():
                ds_save.to_netcdf(f"{sec_path}/{k}.nc")

class BoundedRegion(GriddedRegion):
    def __init__(self, section, grid, **kwargs):
        super().__init__(
            section.name,
            section.lons_c,
            section.lats_c,
            grid,
            **kwargs
        )
        self.children = {}

        def is_slice_in_list(s,l):
            """Check if a list slice (ordered sublist) is contained in list"""
            len_s = len(s) #so we don't recompute length of s on every iteration
            return any(s == l[i:len_s+i] for i in range(len(l) - len_s+1))
            
        for child_name, child in section.children.items():
            i, j, lons, lats = sec.grid_section(grid, child.lons_c, child.lats_c)
            
            lonlat_parent = sec.coords_from_lonlat(loop(self.lons_c), loop(self.lats_c))
            lonlat_child = sec.coords_from_lonlat(lons, lats)
            if not(is_slice_in_list(lonlat_child, lonlat_parent)):
                if is_slice_in_list(lonlat_child[::-1], lonlat_parent):
                    child.lons_c = child.lons_c[::-1]
                    child.lats_c = child.lats_c[::-1]
                    i, j = i[::-1], j[::-1]
                    lons, lats = lons[::-1], lats[::-1]
                else:
                    child.lons_c = child.lons_c[::-1]
                    child.lats_c = child.lats_c[::-1]
                    i, j, lons, lats = sec.grid_section(grid, child.lons_c, child.lats_c)
                    lonlat_parent = sec.coords_from_lonlat(loop(self.lons_c), loop(self.lats_c))
                    lonlat_child = sec.coords_from_lonlat(lons, lats)
                    if not(is_slice_in_list(lonlat_child, lonlat_parent)):
                        raise ValueError("Child sections do not match up with parent ones!")  

            coords = sec.coords_from_lonlat(lons, lats)
            child_section = sec.Section(
                child_name,
                coords,
                children={},
                parents={}
            )
            child_section.i_c = i
            child_section.j_c = j
            self.children[child_name] = child_section

def open_gr(path, ds_to_grid):

    ds_grid = xr.open_dataset(f"{path}/grid.nc")
    grid = ds_to_grid(ds_grid)
    ds = xr.open_dataset(f"{path}/region.nc")
    
    name = path.split('/')[-1][:-3].replace('_',' ')
    region = GriddedRegion(
        name,
        ds.lons_c.values,
        ds.lats_c.values,
        grid,
        mask = ds.mask,
        ij = (ds.i_c.values, ds.j_c.values)
    )
    gr_files = [
        f for f in os.listdir(f"{path}/")
        if ('.nc' in f) and (f not in ['grid.nc', 'region.nc'])
    ]
    for file in gr_files:
        v = file.split('.')[0]
        region.save[v] = xr.open_dataset(f"{path}/{file}")

    region.children = {}
    children_path = f"{path}/children/"
    child_paths = [
        f"{children_path}{file}"
        for file in os.listdir(children_path)
    ]
    for child_path in child_paths:
        child_name = child_path.split('/')[-1][:-4].replace('_', ' ')
        ds = xr.open_dataset(f"{child_path}/section.nc")
        
        section = sec.Section(
            child_name,
            sec.coords_from_lonlat(ds.lons_c.values, ds.lats_c.values)
        )
        
        section.save = {}
        for file in [f for f in os.listdir(child_path) if f != 'section.nc']:
            v = file.split('.')[0]
            section.save[v] = xr.open_dataset(f"{child_path}/{file}")

        region.children[child_name] = section
        
    return region