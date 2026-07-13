import geopandas as gpd
from shapely.geometry import Polygon
import regionmask

import sectionate as sec
from sectionate.gridutils import get_facedim, get_geo_corners, coord_dict
import numpy as np
import xarray as xr
import warnings

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
    """Extend a pole-encircling boundary to the enclosed pole, forming a closed
    polygon that encloses everything between the boundary and that pole.

    The pole is chosen from the boundary's mean latitude: a boundary hugging the
    Arctic encloses the North Pole (+90), one hugging the Antarctic the South Pole
    (-90). The equatorial / ambiguous case (mean latitude 0) defaults to the South
    Pole. The boundary winding is first normalized to eastward, then rolled to be
    single-valued in longitude, then closed off with two points at the pole. Both
    windings of the same boundary therefore enclose the same region.
    """
    pole = 90. if np.mean(lats_c) > 0. else -90.
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
    lats = np.append(lats, [lats[max_idx], pole, pole, lats[min_idx]])

    return Polygon(zip(lons, lats))


def get_region_boundary_grid_indices(lons, lats, grid, curve="latitude circle"):
    """Find boundary coordinates and grid indices that approximate a polygon.

    ARGUMENTS
    ---------
    lons : list or np.ndarray of longitudes
    lats : list or np.ndarray of latitudes
    grid : `xgcm.Grid` instance
    curve : str
        Curve followed between consecutive boundary vertices when snapping them
        onto the grid, passed through to `sectionate.grid_section`. Default:
        ``"latitude circle"`` (constant latitude, marching in longitude), so an
        edge between two vertices at the same latitude stays on that latitude
        rather than bowing poleward as the ``"great circle"`` geodesic would.

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
        sec.grid_section(grid, lons, lats, curve=curve)
    )
    lons_uv, lats_uv = sec.uvcoords_from_qindices(grid, i_c, j_c, f_c=f_c)

    return (i_c, j_c, f_c, lons_c, lats_c, lons_uv, lats_uv)


def rasterize_per_tile(grid, per_slice):
    """Apply a 2D lon/lat rasterization to each tile of a grid and stitch the result.

    ``per_slice(center_lon_2d, center_lat_2d)`` must map a pair of 2D tracer-center
    coordinate DataArrays to a 2D boolean array (or DataArray) of the same shape --
    typically a ``regionmask`` call, which accepts only 1D/2D lon/lat. On a
    single-tile grid (``get_facedim(grid) is None``) it is called once on the full
    2D coordinates. On a multi-tile grid (e.g. a lat-lon-cap / cubed-sphere grid,
    whose coordinates carry an extra face dimension) it is called once per face and
    the per-face masks are stitched back together along the face dimension, so
    ``regionmask`` never sees the 3D coordinate array.

    ARGUMENTS
    ---------
    grid [xgcm.Grid] -- ocean model grid
    per_slice [callable] -- ``(lon_2d, lat_2d) -> 2D bool array/DataArray``

    RETURNS
    -------
    mask : xr.DataArray of bool type, over the grid's tracer-center dims
    """
    geo = get_geo_centers(grid)
    center_lon = geo["X"]
    center_lat = geo["Y"]
    facedim = get_facedim(grid)

    if facedim is None or facedim not in center_lat.dims:
        out = per_slice(center_lon, center_lat)
        return xr.DataArray(
            np.asarray(out, dtype=bool), dims=center_lat.dims, coords=center_lat.coords
        )

    arr = np.zeros(center_lat.shape, dtype=bool)
    for f in range(center_lat.sizes[facedim]):
        idx = tuple(f if d == facedim else slice(None) for d in center_lat.dims)
        with warnings.catch_warnings():
            # A region typically covers only some tiles; regionmask warns "No
            # gridpoint belongs to any region" for the empty ones, which is the
            # expected case when stitching per tile.
            warnings.filterwarnings(
                "ignore", message="No gridpoint belongs to any region"
            )
            out_f = per_slice(
                center_lon.isel({facedim: f}), center_lat.isel({facedim: f})
            )
        arr[idx] = np.asarray(out_f, dtype=bool)
    return xr.DataArray(arr, dims=center_lat.dims, coords=center_lat.coords)


def mask_from_grid_boundaries(
    lons_c,
    lats_c,
    grid,
    ):
    """Find the boolean cell mask bounded by a sequence of cell-corner coordinates.

    Builds a shapely Polygon from ``(lons_c, lats_c)``, splits it at the ±180°
    antimeridian into a clean ``[-180, 180]`` Polygon/MultiPolygon, rasterizes
    each sub-polygon onto the grid's tracer-center lon/lat with regionmask
    (``wrap_lon=False``), and ORs the per-piece boolean masks together. The
    multipolygon pieces "stitch trivially by adding the masks", which handles
    antimeridian-crossing and pole-encircling regions uniformly. Rasterization is
    done per tile via `rasterize_per_tile`, so multi-tile (lat-lon-cap /
    cubed-sphere) grids are supported as well as single-tile ones.

    ARGUMENTS
    ---------
    lons_c [list or np.ndarray] -- cell corner longitudes
    lats_c [list or np.ndarray] -- cell corner latitudes
    grid [xgcm.Grid] -- ocean model grid

    RETURNS
    -------
    region_grid_mask : xr.DataArray of bool type, over the tracer-center dims
    """
    lons_c = np.asarray(lons_c, dtype=float)
    lats_c = np.asarray(lats_c, dtype=float)

    # Total signed longitude winding along the boundary, ignoring antimeridian
    # jumps. A magnitude near 360 means the boundary encircles a pole and cannot
    # be drawn as a simple lon/lat polygon; we then extend it to the enclosed pole
    # so the polygon encloses everything on the boundary's enclosed side.
    dlon = np.diff(lons_c)
    delta_lon = np.sum(dlon[np.abs(dlon) < 180.])

    if np.abs(delta_lon) >= 180.:
        polygon = _pole_enclosing_polygon(lons_c, lats_c, delta_lon)
    else:
        polygon = Polygon(zip(normalize_lon(lons_c), lats_c))

    # Split at the ±180 antimeridian into clean [-180, 180] pieces.
    split = split_at_antimeridian(polygon)

    if split.geom_type == "Polygon":
        pieces = [p for p in [split] if not p.is_empty]
    else:
        pieces = [p for p in split.geoms if not p.is_empty]

    crs = "epsg:4326"

    def _rasterize_pieces(center_lon, center_lat):
        # The grid may use any longitude convention; the polygon pieces live in
        # [-180, 180], so normalize the grid longitudes to match (wrap_lon=False).
        clon = normalize_lon(center_lon)
        piece_mask = np.zeros(center_lat.shape, dtype=bool)
        for piece in pieces:
            gdf = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[piece])
            piece_mask |= ~np.isnan(
                regionmask.mask_geopandas(gdf, clon, lat=center_lat, wrap_lon=False)
            ).values
        return piece_mask

    return rasterize_per_tile(grid, _rasterize_pieces)
