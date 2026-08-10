"""Antimeridian / multipolygon shapely helpers.

These are pure-`shapely` (plus `numpy`) utilities for handling polygons that
cross the ±180° antimeridian or encircle a pole. They are used by
``grid_conform.mask_from_grid_boundaries`` to split a (potentially
antimeridian-crossing) region polygon into a clean ``[-180, 180]`` Polygon or
MultiPolygon whose pieces can each be rasterized independently and OR-ed
together.

Ported from ``xcryocouple.geometry`` (https://github.com/hdrake/xcryocouple).
"""

import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely.geometry import box as _box
from shapely.affinity import translate as _translate
from shapely.validation import make_valid as _make_valid

__all__ = [
    "normalize_lon",
    "has_antimeridian_crossing",
    "polygon_to_lon360",
    "split_at_antimeridian",
]


def normalize_lon(lon):
    """Normalize longitude(s) to the ``[-180, 180]`` range.

    Works on scalars, numpy arrays, and xarray DataArrays. Round-trips: applying
    it twice gives the same result as applying it once.

    Parameters
    ----------
    lon : float, np.ndarray, or xr.DataArray
        Longitude(s) in degrees, in any convention (e.g. ``[0, 360]`` or
        ``[-300, 60]``).

    Returns
    -------
    Same type as input, with longitudes in ``[-180, 180)``.
    """
    return ((lon + 180.) % 360.) - 180.


def has_antimeridian_crossing(geom, threshold=180.0):
    """Return True if the polygon exterior has a longitude jump > *threshold* degrees.

    A large jump between consecutive exterior vertices is the signature of a
    ring that wraps across the ±180° antimeridian when stored in ``[-180, 180]``.

    The default ``threshold`` of 180° cleanly separates real seam crossings from
    legitimate steps when vertices are normalized to ``[-180, 180]``: any single
    grid cell spans < 180° of longitude, so legitimate corner-to-corner steps are
    always < 180°, whereas a seam crossing jumps from near ``+180`` to near
    ``-180`` (magnitude ``360 - cell_width > 180``). Pass a smaller threshold for
    densely-sampled shapefile polygons if desired.

    Parameters
    ----------
    geom : shapely.geometry.Polygon
    threshold : float
        Minimum absolute longitude difference (degrees) flagged as a crossing.

    Returns
    -------
    bool
    """
    if geom.geom_type != "Polygon":
        return False
    lons = np.array([c[0] for c in geom.exterior.coords])
    if lons.size < 2:
        return False
    return bool(np.any(np.abs(np.diff(lons)) > threshold))


def polygon_to_lon360(geom):
    """Return a copy of *geom* with all longitudes converted to ``[0, 360]``.

    Applying ``lon % 360`` to every vertex turns an antimeridian-crossing
    polygon (which has discontinuous ~-180 -> ~+180 jumps in ``[-180, 180]``)
    into a topologically continuous polygon in ``[0, 360]`` space.

    Parameters
    ----------
    geom : shapely.geometry.Polygon

    Returns
    -------
    shapely.geometry.Polygon
    """
    def _ring360(ring):
        return [(lon % 360., lat) for lon, lat in ring.coords]

    exterior = _ring360(geom.exterior)
    interiors = [_ring360(r) for r in geom.interiors]
    return Polygon(exterior, interiors)


def _polygon_pieces(geom):
    """Yield the (non-empty) Polygon components of *geom*."""
    if geom.is_empty:
        return
    if geom.geom_type == "Polygon":
        yield geom
    elif hasattr(geom, "geoms"):
        for g in geom.geoms:
            if g.geom_type == "Polygon" and not g.is_empty:
                yield g


def _as_polygon_or_multipolygon(pieces):
    """Collapse a list of Polygons into a single Polygon or a MultiPolygon."""
    pieces = [p for p in pieces if not p.is_empty]
    if len(pieces) == 1:
        return pieces[0]
    return MultiPolygon(pieces)


def split_at_antimeridian(geom):
    """Split a polygon at ±180° into clean ``[-180, 180]`` pieces.

    Antimeridian-crossing polygons stored in EPSG:4326 have vertices that jump
    discontinuously between ~-180° and ~+180°. Shapely would interpret the ring
    as crossing through the interior of the globe, producing a nearly-global
    bounding box. This function instead converts the polygon to ``[0, 360]``
    space (where the antimeridian sits at 360° and the ring stays continuous),
    ``make_valid``s it there, then splits at 360° and translates any
    ``[180, 360]`` portion back into ``[-180, 0]`` so that all output pieces lie
    within ``[-180, 180]``.

    The antimeridian check is done on the *original* geometry, before
    ``make_valid``, because ``make_valid`` would shred a self-intersecting
    antimeridian ring into many small fragments.

    Parameters
    ----------
    geom : shapely.geometry.Polygon

    Returns
    -------
    shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        A single Polygon when the input does not cross the antimeridian (and is
        already valid as one piece); a MultiPolygon when the split produces
        multiple components. All pieces lie within ``[-180, 180]``.
    """
    if geom.geom_type == "Polygon" and has_antimeridian_crossing(geom):
        # Move to [0, 360] for a topologically continuous ring, then clean.
        geom360 = _make_valid(polygon_to_lon360(geom))

        # Handle the (unlikely) case of also straddling the 360° meridian.
        west360 = geom360.intersection(_box(0, -90, 360, 90))
        east360 = geom360.intersection(_box(360, -90, 720, 90))

        pieces = []
        for part in (west360, east360):
            if part.is_empty:
                continue
            # Left half [0, 180] is already in [-180, 180];
            # right half [180, 360] is translated to [-180, 0].
            left = part.intersection(_box(0, -90, 180, 90))
            right = part.intersection(_box(180, -90, 360, 90))
            if not right.is_empty:
                right = _translate(right, xoff=-360.0)
            pieces.extend(_polygon_pieces(left))
            pieces.extend(_polygon_pieces(right))

        if pieces:
            return _as_polygon_or_multipolygon(pieces)
        return geom

    # Not an antimeridian crossing: just clean it up.
    geom = _make_valid(geom)
    pieces = list(_polygon_pieces(geom))
    if pieces:
        return _as_polygon_or_multipolygon(pieces)
    return geom
