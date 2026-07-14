import numpy as np
import xarray as xr

import sectionate as sec
from sectionate import is_section_counterclockwise
from sectionate.gridutils import get_geo_corners
from .utilities import *
from .grid_conform import (
    get_region_boundary_grid_indices,
    mask_from_grid_boundaries,
    _normalize_grid_section,
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
        curve="great circle",
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
        curve : str
            Curve followed between consecutive boundary vertices when snapping the
            polygon onto the grid. Default: ``"great circle"`` (the geodesic), which
            traces arbitrary edges -- including diagonal ones -- correctly. Pass
            ``"latitude circle"`` to instead march along constant latitude, so a box
            edge between two same-latitude vertices stays on that latitude rather than
            bowing poleward (only sensible for axis-aligned boxes; it snaps a diagonal
            edge onto a latitude line). Only used on the boundary-defined construction
            path (ignored when `ij` is given).

        RETURNS
        -------
        GriddedRegion instance

        Examples
        --------
        >>> grid = xgcm.Grid(...) # TO DO: minimal example grid
        >>> lons, lats = np.array([-80., -66., -65.]), np.array([ 26.,  18.,  32.])
        >>> region = reg.Region('Bermuda Triangle', lons, lats, grid)
        """

        try:
            get_geo_corners(grid)
        except ValueError as e:
            raise ValueError(
                "grid._ds must contain two-dimensional cell-corner (vorticity) "
                'coordinates whose names contain "lon" and "lat".'
            ) from e

        self.grid = grid
        self.save = {}

        if len(lons)>=3 and len(lats)>=3 and ij is None:
            self.initiate_from_boundary(
                lons,
                lats,
                mask=mask,
                positive_in=positive_in,
                curve=curve
            )
        elif ij is None:
            raise NameError("Must provide lons and lats as lists or arrays\
            to define the region.")
        else:
            self.lons_c = lons
            self.lats_c = lats
            self.i_c = ij[0]
            self.j_c = ij[1]
            self.f_c = ij[2] if len(ij) > 2 else None
            if mask is None:
                self.mask = mask_from_grid_boundaries(
                    self.lons_c,
                    self.lats_c,
                    self.grid,
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
        mask=None,
        curve="great circle"
        ):
        """
        TO DO
        """

        (self.i_c, self.j_c, self.f_c, self.lons_c, self.lats_c,
         self.lons_uv, self.lats_uv) = (
            get_region_boundary_grid_indices(
                lons.copy(),
                lats.copy(),
                self.grid,
                curve=curve
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
          - a `region.nc` file that contains the region's tracer cell mask
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
        >>> gridded_region.to_gr('../data/')
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

        # Write boundary information to NetCDF file. The corner indices
        # (i_c, j_c[, f_c]) and the boundary coordinates (lons_c, lats_c) may
        # have different lengths (the indices are closed loops; the coordinates
        # are open), so they live on separate dimensions.
        ds = xr.Dataset()
        ds['lons_c'] = xr.DataArray(np.asarray(self.lons_c), dims=('vertex',))
        ds['lats_c'] = xr.DataArray(np.asarray(self.lats_c), dims=('vertex',))
        ds['i_c'] = xr.DataArray(np.asarray(self.i_c), dims=('corner',))
        ds['j_c'] = xr.DataArray(np.asarray(self.j_c), dims=('corner',))
        if getattr(self, 'f_c', None) is not None:
            ds['f_c'] = xr.DataArray(np.asarray(self.f_c), dims=('corner',))
        for v in ['lons_uv', 'lats_uv']:
            if getattr(self, v, None) is not None:
                ds[v] = xr.DataArray(np.asarray(getattr(self, v)), dims=('face',))
        ds['mask'] = self.mask
        ds.to_netcdf(f"{gr_path}/region.nc")

        for (k,v) in self.save.items():
            v.to_netcdf(f"{gr_path}/{k}.nc")

        # Write boundary information for each child section to NetCDF file
        children = getattr(self, 'children', {})
        if children:
            child_path = f"{gr_path}/children/"
            Path(child_path).mkdir(parents=True, exist_ok=True)
            for child in children.values():
                child_name = child.name.replace(' ','_')
                sec_path = f"{gr_path}/children/{child_name}.sec"
                Path(sec_path).mkdir(parents=True, exist_ok=True)

                ds = xr.Dataset()
                ds['lons_c'] = xr.DataArray(np.asarray(child.lons_c), dims=('vertex',))
                ds['lats_c'] = xr.DataArray(np.asarray(child.lats_c), dims=('vertex',))
                ds['i_c'] = xr.DataArray(np.asarray(child.i_c), dims=('corner',))
                ds['j_c'] = xr.DataArray(np.asarray(child.j_c), dims=('corner',))
                if getattr(child, 'f_c', None) is not None:
                    ds['f_c'] = xr.DataArray(np.asarray(child.f_c), dims=('corner',))
                ds.to_netcdf(f"{sec_path}/section.nc")
                for k, ds_save in getattr(child, 'save', {}).items():
                    ds_save.to_netcdf(f"{sec_path}/{k}.nc")

class BoundedRegion(GriddedRegion):
    def __init__(self, section, grid, curve="great circle", **kwargs):
        super().__init__(
            section.name,
            section.lons_c,
            section.lats_c,
            grid,
            curve=curve,
            **kwargs
        )
        self.children = {}

        def slice_indices_in_list(s, l):
            """
            Return a list of indices in `l` where elements of slice `s` appear in order.
            Returns list of indices if found, else None.
            """
            len_s = len(s)
            for i in range(len(l) - len_s + 1):
                if s == l[i:i+len_s]:
                    return np.array(list(range(i, i+len_s)), dtype=int)
            return None

        # Create the parent gridded section
        parent_section_gridded = sec.GriddedSection(
            sec.Section(
                section.name,
                (self.lons_c, self.lats_c),
            ),
            grid,
            i_c = self.i_c,
            j_c = self.j_c,
            f_c = self.f_c,
        )
        parent_lons_uv, parent_lats_uv = sec.uvcoords_from_qindices(
            grid,
            parent_section_gridded.i_c,
            parent_section_gridded.j_c,
            f_c = parent_section_gridded.f_c,
        )
        parent_coords_uv = sec.coords_from_lonlat(parent_lons_uv, parent_lats_uv)
            
        for child_name, child in section.children.items():
            i_c, j_c, f_c, lons_c, lats_c = _normalize_grid_section(
                sec.grid_section(grid, child.lons_c, child.lats_c, curve=curve)
            )

            child_coords = sec.coords_from_lonlat(lons_c, lats_c)

            # Find the indices in the parent's corner sections that correspond this child
            parent_idx_c = slice_indices_in_list(child_coords, parent_section_gridded.coords)
            if parent_idx_c is not None:
                pass  # child orientation matches parent
            else:
                parent_idx_c = slice_indices_in_list(child_coords[::-1], parent_section_gridded.coords)
                if parent_idx_c is not None:
                    # reversed slice is in the parent list
                    child.lons_c = child.lons_c[::-1]
                    child.lats_c = child.lats_c[::-1]
                    # recompute the child sections using the correct orientation
                    i_c, j_c, f_c, lons_c, lats_c = _normalize_grid_section(
                        sec.grid_section(grid, child.lons_c, child.lats_c, curve=curve)
                    )
                else:
                    raise ValueError("Child corner sections do not match up with parent ones!")

            # Find the indices in the parent's velocity sections that correspond this child
            lons_uv, lats_uv = sec.uvcoords_from_qindices(grid, i_c, j_c, f_c=f_c)
            child_coords_uv = sec.coords_from_lonlat(lons_uv, lats_uv)
            parent_idx_uv = slice_indices_in_list(child_coords_uv, parent_coords_uv)
            if parent_idx_uv is None:
                raise ValueError("Child velocity sections do not match up with parent ones!")

            child_coords = sec.coords_from_lonlat(lons_c, lats_c)
            child_section = sec.Section(
                child_name,
                child_coords,
                children={},
                parent=parent_section_gridded
            )
            child_section_gridded = sec.GriddedSection(
                child_section,
                grid,
                i_c = i_c,
                j_c = j_c,
                f_c = f_c,
            )
            child_section_gridded.parent_idx_c = parent_idx_c
            child_section_gridded.parent_idx_uv = parent_idx_uv
            self.children[child_name] = child_section_gridded

class MaskRegion:
    """One topology-aware connected component of a cell mask on a C-grid model.

    Unlike `GriddedRegion` -- a single closed polygon boundary that is *handed* a
    mask -- a `MaskRegion` owns exactly the cells of one connected component
    (`.mask`, derived from the labeling, hence unambiguous) together with the full
    set of grid-conforming boundary loops that enclose *those* cells
    (`.boundaries`, each a `sectionate.GriddedSection`). A component that wraps a
    seam or contains holes simply has more than one boundary loop; integrating a
    flux over every loop reproduces the flux convergence over `.mask` exactly (the
    discrete divergence theorem), whatever the loop count.

    PARAMETERS
    ----------
    name : str
    grid : `xgcm.Grid` instance
    mask : `xr.DataArray` of bool -- this component's own cells
    boundaries : list of `sectionate.GriddedSection` -- the loops enclosing `mask`
    """
    def __init__(self, name, grid, mask, boundaries):
        self.name = name
        self.grid = grid
        self.mask = mask
        self.boundaries = boundaries
        self.save = {}

    def __repr__(self):
        n = len(self.boundaries)
        return (f"{str(type(self))[8:-2]}('{self.name}', "
                f"{n} boundar{'y' if n == 1 else 'ies'})")

    def to_gr(self, path):
        """Save the MaskRegion as a `.gr` directory.

        Unlike `GriddedRegion.to_gr` (a single boundary loop stored in
        ``region.nc``), a `MaskRegion` has a *list* of boundary loops, so each is
        written as its own sub-section under ``boundaries/``. The directory is
        tagged ``kind="MaskRegion"`` on ``region.nc`` so `open_gr` knows to rebuild
        a `MaskRegion` rather than a `GriddedRegion`.

        Layout::

            <name>.gr/
              grid.nc                            # grid coords (or symlink to ../grid.nc)
              region.nc                          # the component's boolean mask
              boundaries/loop_<k>.sec/section.nc # each loop's lons_c/lats_c/i_c/j_c[/f_c]

        Arguments
        ---------
        path [str] -- directory to write ``[MaskRegion.name].gr`` into.
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

        ds = xr.Dataset(attrs={"kind": "MaskRegion"})
        ds['mask'] = self.mask
        ds.to_netcdf(f"{gr_path}/region.nc")

        for (k, v) in self.save.items():
            v.to_netcdf(f"{gr_path}/{k}.nc")

        bnd_path = f"{gr_path}/boundaries/"
        Path(bnd_path).mkdir(parents=True, exist_ok=True)
        for k, loop in enumerate(self.boundaries):
            sec_path = f"{bnd_path}/loop_{k}.sec"
            Path(sec_path).mkdir(parents=True, exist_ok=True)
            dsb = xr.Dataset()
            dsb['lons_c'] = xr.DataArray(np.asarray(loop.lons_c), dims=('vertex',))
            dsb['lats_c'] = xr.DataArray(np.asarray(loop.lats_c), dims=('vertex',))
            dsb['i_c'] = xr.DataArray(np.asarray(loop.i_c), dims=('corner',))
            dsb['j_c'] = xr.DataArray(np.asarray(loop.j_c), dims=('corner',))
            if getattr(loop, 'f_c', None) is not None:
                dsb['f_c'] = xr.DataArray(np.asarray(loop.f_c), dims=('corner',))
            dsb.to_netcdf(f"{sec_path}/section.nc")


def _open_mask_region_gr(path, name, grid, ds):
    """Reconstruct a `MaskRegion` from a ``.gr`` directory written by
    `MaskRegion.to_gr` (its mask plus one gridded-section loop per ``boundaries/``
    sub-directory, each carrying its stored i_c/j_c[/f_c])."""
    mask = ds['mask']
    bnd_path = f"{path}/boundaries/"
    loop_dirs = sorted(
        [d for d in os.listdir(bnd_path) if d.endswith('.sec')],
        key=lambda d: int(d[len('loop_'):-len('.sec')]),
    )
    boundaries = []
    for d in loop_dirs:
        dsb = xr.open_dataset(f"{bnd_path}/{d}/section.nc")
        f_c = dsb.f_c.values if 'f_c' in dsb else None
        boundaries.append(sec.GriddedSection(
            sec.Section(d[:-4], sec.coords_from_lonlat(dsb.lons_c.values, dsb.lats_c.values)),
            grid, i_c=dsb.i_c.values, j_c=dsb.j_c.values, f_c=f_c,
        ))
    region = MaskRegion(name, grid, mask, boundaries)
    for file in [f for f in os.listdir(path)
                 if f.endswith('.nc') and f not in ('grid.nc', 'region.nc')]:
        region.save[file.split('.')[0]] = xr.open_dataset(f"{path}/{file}")
    return region


def open_gr(path, ds_to_grid):

    ds_grid = xr.open_dataset(f"{path}/grid.nc")
    grid = ds_to_grid(ds_grid)
    ds = xr.open_dataset(f"{path}/region.nc")

    name = path.split('/')[-1][:-3].replace('_',' ')

    # A MaskRegion `.gr` stores only the mask in region.nc (its loops live under
    # boundaries/); a GriddedRegion `.gr` stores a single loop's i_c/lons_c inline.
    if ds.attrs.get('kind') == 'MaskRegion' or os.path.isdir(f"{path}/boundaries"):
        return _open_mask_region_gr(path, name, grid, ds)

    f_c = ds.f_c.values if 'f_c' in ds else None
    region = GriddedRegion(
        name,
        ds.lons_c.values,
        ds.lats_c.values,
        grid,
        mask = ds.mask,
        ij = (ds.i_c.values, ds.j_c.values, f_c)
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
    if not os.path.isdir(children_path):
        return region
    child_paths = [
        f"{children_path}{file}"
        for file in os.listdir(children_path)
    ]
    for child_path in child_paths:
        child_name = child_path.split('/')[-1][:-4].replace('_', ' ')
        ds = xr.open_dataset(f"{child_path}/section.nc")

        # reconstruct the child as a gridded section carrying its stored corner
        # indices (i_c/j_c/f_c), mirroring how `BoundedRegion` builds children --
        # rather than discarding them and rebuilding a bare `sec.Section` from coords.
        child_f_c = ds.f_c.values if 'f_c' in ds else None
        section = sec.GriddedSection(
            sec.Section(
                child_name,
                sec.coords_from_lonlat(ds.lons_c.values, ds.lats_c.values),
            ),
            grid,
            i_c=ds.i_c.values,
            j_c=ds.j_c.values,
            f_c=child_f_c,
        )

        section.save = {}
        for file in [f for f in os.listdir(child_path) if f != 'section.nc']:
            v = file.split('.')[0]
            section.save[v] = xr.open_dataset(f"{child_path}/{file}")

        region.children[child_name] = section
        
    return region