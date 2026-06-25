import geopandas as gpd
from shapely.geometry import Polygon
import regionmask

import sectionate as sec
from sectionate.gridutils import get_facedim, get_geo_corners, coord_dict
import numpy as np

from .utilities import *
from .geometry import split_at_antimeridian, normalize_lon


def _normalize_grid_section(result):
    """Normalize `sec.grid_section`'s 4-or-5-tuple return to always include `f_c`.

    `sec.grid_section` returns ``(i_c, j_c, lons_c, lats_c)`` for single-tile grids
    and ``(i_c, j_c, f_c, lons_c, lats_c)`` for multi-tile grids. Downstream code
    always wants ``(i_c, j_c, f_c, lons_c, lats_c)`` with ``f_c`` set to None for
    single-tile grids.
    """
    if len(result) == 5:
        i_c, j_c, f_c, lons_c, lats_c = result
    elif len(result) == 4:
        i_c, j_c, lons_c, lats_c = result
        f_c = None
    else:
        raise ValueError(
            f"Unexpected `grid_section` return of length {len(result)}; "
            "expected a 4-tuple (single-tile) or 5-tuple (multi-tile)."
        )
    return i_c, j_c, f_c, lons_c, lats_c


def get_geo_centers(grid):
    """Find the tracer-center longitude and latitude coordinate DataArrays.

    Analogous to `sectionate.gridutils.get_geo_corners`, but returns the
    coordinates at the cell *center* (tracer) position rather than the corner
    (vorticity) position. The center dimension names are discovered via
    `coord_dict(grid)`, and the matching ``"lon"``/``"lat"`` coordinate
    variables are selected by substring -- so no coordinate names are
    hard-coded.

    Parameters
    ----------
    grid : xgcm.Grid

    Returns
    -------
    dict
        ``{"X": <center-lon DataArray>, "Y": <center-lat DataArray>}``.
    """
    cdict = coord_dict(grid)
    Xdim = cdict["X"]["center"]
    Ydim = cdict["Y"]["center"]

    coords = grid._ds.coords
    geo = {}
    for axis, geoc in zip(["X", "Y"], ["lon", "lat"]):
        matches = [
            coords[c] for c in coords
            if (geoc in c.lower())
            and (Xdim in coords[c].dims)
            and (Ydim in coords[c].dims)
        ]
        if len(matches) == 0:
            raise ValueError(
                'grid._ds must contain two-dimensional ("X", "Y") tracer-center '
                'coordinates including the strings "lon" and "lat".'
            )
        geo[axis] = matches[0]
    return geo


def get_region_boundary_grid_indices(lons, lats, grid):
    """Find boundary coordinates and grid indices that approximate a polygon.

    ARGUMENTS
    ---------
    lons : list or np.ndarray of longitudes
    lats : list or np.ndarray of latitudes
    grid : `xgcm.Grid` instance

    RETURNS
    -------
    (i_c, j_c, f_c, lons_c, lats_c, lons_uv, lats_uv)

    i_c : "X"-axis grid indices of corner points
    j_c : "Y"-axis grid indices of corner points
    f_c : face/tile indices of corner points (None for single-tile grids)
    lons_c : longitudes of corner points
    lats_c : latitudes of corner points
    lons_uv : longitudes of (u,v) velocity faces
    lats_uv : latitudes of (u,v) velocity faces
    """
    if (lons[0], lats[0]) != (lons[-1], lats[-1]):
        lons, lats = loop(lons), loop(lats)

    i_c, j_c, f_c, lons_c, lats_c = _normalize_grid_section(
        sec.grid_section(grid, lons, lats)
    )
    lons_uv, lats_uv = sec.uvcoords_from_qindices(grid, i_c, j_c, f_c=f_c)

    return (i_c, j_c, f_c, lons_c, lats_c, lons_uv, lats_uv)


def mask_from_grid_boundaries(
    lons_c,
    lats_c,
    grid,
    along_boundary=False,
    ):
    """Find the boolean cell mask bounded by a sequence of cell-corner coordinates.

    Builds a shapely Polygon from ``(lons_c, lats_c)``, splits it at the ±180°
    antimeridian into a clean ``[-180, 180]`` Polygon/MultiPolygon, rasterizes
    each sub-polygon onto the grid's tracer-center lon/lat with regionmask
    (``wrap_lon=False``), and ORs the per-piece boolean masks together. The
    multipolygon pieces "stitch trivially by adding the masks", which handles
    antimeridian-crossing and pole-encircling regions uniformly.

    ARGUMENTS
    ---------
    lons_c [list or np.ndarray] -- cell corner longitudes
    lats_c [list or np.ndarray] -- cell corner latitudes
    grid [xgcm.Grid] -- ocean model grid
    along_boundary [bool] -- if True, treat (lons_c, lats_c) as an already-closed
        boundary polygon (no extra closure logic). Currently both paths use the
        same split+OR rasterization; the flag is retained for API compatibility.

    RETURNS
    -------
    region_grid_mask : xr.DataArray of bool type, over the tracer-center dims
    """
    geo = get_geo_centers(grid)
    center_lon = normalize_lon(geo["X"])
    center_lat = geo["Y"]

    # Build a polygon from the (normalized) corner coordinates and split it at
    # the antimeridian into clean [-180, 180] pieces.
    lons_n = normalize_lon(np.asarray(lons_c, dtype=float))
    lats_n = np.asarray(lats_c, dtype=float)
    polygon = Polygon(zip(lons_n, lats_n))
    split = split_at_antimeridian(polygon)

    if split.geom_type == "Polygon":
        pieces = [split]
    else:
        pieces = list(split.geoms)

    crs = "epsg:4326"
    region_grid_mask = None
    for piece in pieces:
        if piece.is_empty:
            continue
        gdf = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[piece])
        piece_mask = ~np.isnan(
            regionmask.mask_geopandas(
                gdf,
                center_lon,
                lat=center_lat,
                wrap_lon=False,
            )
        )
        if region_grid_mask is None:
            region_grid_mask = piece_mask
        else:
            region_grid_mask = region_grid_mask | piece_mask

    return region_grid_mask
