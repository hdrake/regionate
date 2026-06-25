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


def _pole_enclosing_polygon(lons_c, lats_c, delta_lon):
    """Extend a pole-encircling boundary down to the South Pole, forming a closed
    polygon that encloses everything south of the boundary.

    The boundary winding is first normalized to eastward (sectionate's
    stereographic-plane orientation convention relative to the South Pole), then
    rolled to be single-valued in longitude, then closed off with two points at
    latitude -90. Both windings of the same boundary therefore enclose the same
    region.
    """
    s = np.sign(delta_lon).astype(int)
    if s == -1:
        lons_c = lons_c[::-1]
        lats_c = lats_c[::-1]
        s = 1

    min_idx = np.argmin(lons_c)
    lons = np.roll(lons_c, -min_idx)
    lats = np.roll(lats_c, -min_idx)

    lons = np.append(lon_mod(lons[-1], lons[0]), lons)
    lats = np.append(lats[-1], lats)

    diffs = s * (lons[np.newaxis, :] - lons[:, np.newaxis])
    diffs[np.tril_indices(lons.size)] *= -1
    single_valued = ~np.any(diffs < 0, axis=1)

    roll_idx = np.argmax(single_valued[::s])
    lons = np.roll(lons[::s], -roll_idx)[::s]
    lats = np.roll(lats[::s], -roll_idx)[::s]

    min_idx = np.argmin(lons)
    max_idx = np.argmax(lons)
    lons = np.append(
        lons, [lons[max_idx] + 10, lons[max_idx] + 10,
               lons[min_idx] - 10, lons[min_idx] - 10]
    )
    lats = np.append(lats, [lats[max_idx], -90, -90, lats[min_idx]])

    return Polygon(zip(lons, lats))


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

    lons_c = np.asarray(lons_c, dtype=float)
    lats_c = np.asarray(lats_c, dtype=float)

    # Total signed longitude winding along the boundary, ignoring antimeridian
    # jumps. A magnitude near 360 means the boundary encircles a pole and cannot
    # be drawn as a simple lon/lat polygon; we then extend it to the South Pole
    # (sectionate's stereographic-plane orientation convention) so the polygon
    # encloses everything on the boundary's enclosed side.
    dlon = np.diff(lons_c)
    delta_lon = np.sum(dlon[np.abs(dlon) < 180.])

    if (not along_boundary) and (np.abs(delta_lon) >= 180.):
        polygon = _pole_enclosing_polygon(lons_c, lats_c, delta_lon)
    else:
        polygon = Polygon(zip(normalize_lon(lons_c), lats_c))

    # Split at the ±180 antimeridian into clean [-180, 180] pieces.
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
