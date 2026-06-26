import contourpy
import numpy as np
import xarray as xr

from .utilities import loop
from sectionate.gridutils import (
    get_facedim,
    get_geo_corners,
    coord_dict,
    check_symmetric,
    build_neighbor_maps,
)
from sectionate.section import distance_on_unit_sphere, COINCIDENT_TOLERANCE_M


def grid_boundaries_from_mask(grid, mask):
    """Find the cell-corner boundaries that enclose a boolean cell `mask`.

    For single-tile grids (no `face_connections`) this traces the mask boundary
    with `contourpy`, exactly as before. For multi-tile grids (lat-lon-cap,
    cubed-sphere, ...) the boundary is traced face-by-face and stitched across
    tile seams using the grid's `xgcm` topology, so a region that spans several
    faces yields a single closed boundary loop.

    Returns lists with a common length equal to the number of discrete contours
    that bound the mask.

    ARGUMENTS
    ---------
    grid : `xgcm.Grid` instance
    mask : `xr.DataArray` instance of type `bool`

    RETURNS
    -------
    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list

    `f_c_list` holds the per-corner face index for multi-tile grids; its entries
    are `None` for single-tile grids.
    """
    if get_facedim(grid) is None:
        return _single_tile_boundaries_from_mask(grid, mask)
    else:
        return _multitile_boundaries_from_mask(grid, mask)


def _contour_to_corner_indices(c, symmetric):
    """Map one `contourpy` polyline (closed, first==last) to closed cell-corner
    index arrays `(i_c, j_c)`. This is the original single-tile remap, factored
    out so it can be shared."""
    i_c, j_c = c[:-1, 0], c[:-1, 1]

    i_c_new, j_c_new = i_c.copy(), j_c.copy()

    i_inc = np.roll(i_c, -1) - i_c
    j_inc = np.roll(j_c, -1) - j_c

    i_c_new[(i_c % 1) == 0.0] = (i_c - (i_inc < 0))[(i_c % 1) == 0.0] + symmetric
    j_c_new[(j_c % 1) == 0.0] = (j_c - (j_inc < 0))[(j_c % 1) == 0.0] + symmetric
    i_c_new[(i_c % 1) == 0.5] = np.floor(i_c[(i_c % 1) == 0.5]) + symmetric
    j_c_new[(j_c % 1) == 0.5] = np.floor(j_c[(j_c % 1) == 0.5]) + symmetric

    return loop(i_c_new).astype(np.int64), loop(j_c_new).astype(np.int64)


def _single_tile_boundaries_from_mask(grid, mask):
    """Trace mask boundaries on a single-tile grid with `contourpy`."""
    symmetric = check_symmetric(grid)
    cdict = coord_dict(grid)
    Xc, Yc = cdict["X"]["center"], cdict["Y"]["center"]
    Xq, Yq = cdict["X"]["corner"], cdict["Y"]["corner"]
    geo = get_geo_corners(grid)
    lon_c, lat_c = geo["X"], geo["Y"]

    m = mask.transpose(Yc, Xc).values

    contours = contourpy.contour_generator(
        np.arange(-1, m.shape[1] + 1),
        np.arange(-1, m.shape[0] + 1),
        np.pad(m, np.array([1, 1])),
    ).create_contour(0.5)

    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list = [], [], [], [], []
    for c in contours:
        i_c_new, j_c_new = _contour_to_corner_indices(c, symmetric)

        i_c_list.append(i_c_new)
        j_c_list.append(j_c_new)
        f_c_list.append(None)

        idx = {
            Xq: xr.DataArray(i_c_new, dims=("pt",)),
            Yq: xr.DataArray(j_c_new, dims=("pt",)),
        }
        lons_c_list.append(lon_c.isel(idx).values[:-1])
        lats_c_list.append(lat_c.isel(idx).values[:-1])

    return i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list


def _pad_mask(grid, mask):
    """Pad the boolean `mask` (on tracer-center points) by one cell in each
    direction using the grid's own `xgcm` topology, so the halo holds each edge
    cell's real neighbour across periodic wraps, tile seams and folds. Genuine
    walls (fill/extend edges) are filled with NaN. Mirrors the padding used by
    `sectionate.gridutils.build_neighbor_maps`, but for center-point data."""
    boundary = {ax: grid.axes[ax].boundary for ax in grid.axes}
    boundary_width = {ax: (1, 1) for ax in grid.axes}
    a = mask.astype(float)
    if hasattr(grid, "pad"):
        return grid.pad(a, boundary_width=boundary_width, boundary=boundary,
                        fill_value=np.nan)
    from xgcm.padding import pad as _module_pad
    return _module_pad(a, grid, boundary_width, boundary=boundary,
                       fill_value=np.nan)


# Each in-mask cell (f, j, i) contributes a boundary face for any of its four
# neighbours that is not in the mask. A boundary face is recorded as a directed
# corner edge oriented so that the in-mask cell lies to its LEFT (counter-
# clockwise winding). Corner indices are expressed relative to the cell using
# the symmetric/non-symmetric offset `o` (o=1 for 'outer' symmetric grids,
# o=0 for 'right' non-symmetric grids): the cell's SW corner is (i-1+o, j-1+o).
_DIRECTED_FACE_EDGES = {
    # neighbour direction: (start_corner_offset, end_corner_offset)
    # corner offsets are (di, dj) added to (i-1+o, j-1+o)
    "right": ((1, 0), (1, 1)),  # east  face: SE -> NE  (+j)
    "left":  ((0, 1), (0, 0)),  # west  face: NW -> SW  (-j)
    "up":    ((1, 1), (0, 1)),  # north face: NE -> NW  (-i)
    "down":  ((0, 0), (1, 0)),  # south face: SW -> SE  (+i)
}


def _multitile_boundaries_from_mask(grid, mask):
    """Trace and stitch mask boundaries across tile seams on a multi-tile grid."""
    facedim = get_facedim(grid)
    symmetric = check_symmetric(grid)
    o = 1 if symmetric else 0

    cdict = coord_dict(grid)
    Xc, Yc = cdict["X"]["center"], cdict["Y"]["center"]
    Xq, Yq = cdict["X"]["corner"], cdict["Y"]["corner"]
    geo = get_geo_corners(grid)
    lon_c, lat_c = geo["X"].transpose(facedim, Yq, Xq), geo["Y"].transpose(facedim, Yq, Xq)

    mask = mask.transpose(facedim, Yc, Xc)
    nf = mask.sizes[facedim]

    Mpad = _pad_mask(grid, mask).transpose(facedim, ..., Yc, Xc).values
    m = mask.values.astype(bool)
    lon_v, lat_v = lon_c.values, lat_c.values

    # Neighbour value of each interior cell across each face (NaN == wall).
    nbr = {
        "right": Mpad[:, 1:-1, 2:],
        "left":  Mpad[:, 1:-1, :-2],
        "up":    Mpad[:, 2:, 1:-1],
        "down":  Mpad[:, :-2, 1:-1],
    }

    def corner(f, ci, cj):
        return (int(f), int(ci), int(cj), float(lon_v[f, cj, ci]), float(lat_v[f, cj, ci]))

    edges = []  # list of (start_corner, end_corner)
    fcells, jcells, icells = np.where(m)
    for f, j, i in zip(fcells, jcells, icells):
        base_i, base_j = i - 1 + o, j - 1 + o
        for d, ((sdi, sdj), (edi, edj)) in _DIRECTED_FACE_EDGES.items():
            neighbour = nbr[d][f, j, i]
            if neighbour == 1.0:  # in-mask neighbour (incl. across a seam): internal
                continue
            start = corner(f, base_i + sdi, base_j + sdj)
            end = corner(f, base_i + edi, base_j + edj)
            edges.append((start, end))

    # Topology-aware corner adjacency (handles rotated/reversed tile seams, where
    # the two tiles' corners are offset and do not coincide in lon/lat).
    maps = build_neighbor_maps(grid, geo)

    def neighbor_corners(f, ci, cj):
        nbrs = set()
        for d in maps:
            fmap, jmap, imap = maps[d]
            nf_ = int(fmap[f, cj, ci]) if fmap is not None else f
            nbrs.add((nf_, int(imap[f, cj, ci]), int(jmap[f, cj, ci])))
        return nbrs

    return _stitch_edges_into_loops(edges, neighbor_corners)


def _stitch_edges_into_loops(edges, neighbor_corners=None):
    """Stitch directed corner edges into ordered closed boundary loops.

    The end corner of each edge is matched to the start corner of the next edge
    in three escalating ways: (1) the exact same grid corner (within a tile);
    (2) a physically coincident corner (a shared seam corner of a non-rotated
    tile connection); and (3) a grid-adjacent corner via the topology neighbour
    maps (a rotated/reversed seam crossing, where the two tiles' corners are
    offset by a cell and never coincide in lon/lat). `neighbor_corners(f, i, j)`
    returns the set of grid-neighbour corners of a corner; pass it for multi-tile
    grids."""
    if not edges:
        return [], [], [], [], []

    # Index edges by their exact start corner (face, i, j).
    start_by_corner = {}
    for k, (s, e) in enumerate(edges):
        start_by_corner.setdefault((s[0], s[1], s[2]), []).append(k)

    used = [False] * len(edges)

    def find_next(end):
        # 1. an edge starting at the exact same corner (interior of a tile)
        for k in start_by_corner.get((end[0], end[1], end[2]), []):
            if not used[k]:
                return k
        # 2. an edge starting at a physically coincident corner (non-rotated seam)
        for k, (s, _e) in enumerate(edges):
            if (not used[k] and
                    distance_on_unit_sphere(end[3], end[4], s[3], s[4]) < COINCIDENT_TOLERANCE_M):
                return k
        # 3. an edge starting at a grid-adjacent corner (rotated/reversed seam)
        if neighbor_corners is not None:
            for nb in neighbor_corners(end[0], end[1], end[2]):
                for k in start_by_corner.get(nb, []):
                    if not used[k]:
                        return k
        return None

    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list = [], [], [], [], []

    for start_k in range(len(edges)):
        if used[start_k]:
            continue
        # Walk a single closed loop starting from this unused edge.
        loop_edges = []
        k = start_k
        while k is not None and not used[k]:
            used[k] = True
            loop_edges.append(k)
            k = find_next(edges[k][1])

        seq = _loop_corner_sequence([edges[k] for k in loop_edges])
        i_c = np.array([c[1] for c in seq], dtype=np.int64)
        j_c = np.array([c[2] for c in seq], dtype=np.int64)
        f_c = np.array([c[0] for c in seq], dtype=np.int64)
        lons = np.array([c[3] for c in seq])
        lats = np.array([c[4] for c in seq])

        i_c_list.append(loop(i_c))
        j_c_list.append(loop(j_c))
        f_c_list.append(loop(f_c))
        lons_c_list.append(lons)
        lats_c_list.append(lats)

    return i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list


def _loop_corner_sequence(loop_edges):
    """Build the open corner sequence for one closed loop of directed edges.

    Each edge's start corner is emitted. Where an edge ends on one face and the
    next edge starts on another (a seam crossing), the far-face end corner is
    emitted too, producing the consecutive duplicate physical points that
    `sectionate.uvindices_from_qindices` consumes (it drops the zero-length
    seam faces)."""
    seq = []
    m = len(loop_edges)
    for k in range(m):
        s, e = loop_edges[k]
        ns = loop_edges[(k + 1) % m][0]  # next edge's start corner
        seq.append(s)
        if (e[0], e[1], e[2]) != (ns[0], ns[1], ns[2]):
            seq.append(e)
    return seq
