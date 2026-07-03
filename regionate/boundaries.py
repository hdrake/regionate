import contourpy
import numpy as np
import xarray as xr
from xgcm.padding import pad
try:
    # north-fold boundary detector; only present in xgcm with north-fold support (xgcm#711)
    from xgcm.padding import _is_fold_boundary
except ImportError:  # pragma: no cover - fall back so `import regionate` works on any xgcm
    from collections.abc import Mapping

    def _is_fold_boundary(boundary):
        return isinstance(boundary, Mapping) and "fold" in boundary

from .utilities import loop
from sectionate.gridutils import (
    get_facedim,
    get_geo_corners,
    coord_dict,
    corner_offset,
    build_neighbor_maps,
    NEIGHBOR_DIRECTIONS,
)


def grid_boundaries_from_mask(grid, mask):
    """Find the cell-corner boundaries that enclose a boolean cell `mask`.

    The tracer follows the grid's own topology, so it stitches a region across every
    seam its `xgcm.Grid` declares -- periodic axes, the bipolar/tripolar north
    **fold** (`boundary={..., "Y": {"fold": ...}}`), and multi-tile `face_connections`
    (lat-lon-cap / cubed-sphere). A region that wraps a seam therefore yields a single
    seam-consistent boundary rather than pieces split at the seam with spurious seam
    faces.

    All cases share one front-end (`_trace_and_drop`): `contourpy` traces the mask per
    face, and boundary segments lying on a seam between two in-mask cells are dropped
    as interior (using the grid's own topology-aware halo, `_pad_center`). The surviving
    arcs are then stitched into closed loops by one of two back-ends:

    - single-tile (periodic and/or walled, incl. the fold) -- stitched by physical
      coincidence at the seam; face index `f_c` is ``None``;
    - multi-tile (`face_connections`) -- stitched by cell-set + grid topology, returning
      a per-corner face index `f_c` (needed for rotated/reversed seams whose two sides
      do not share coordinates).

    In every case the returned loops' corners map cleanly to velocity faces via
    `sectionate.uvindices_from_qindices`, so integrating a flux over the boundary
    reproduces its convergence over the masked cells exactly (a property of the full
    loop *set* -- an annulus, fold, or multi-tile region may return several loops).

    Returns lists with a common length equal to the number of discrete boundary loops.
    `f_c_list` holds the per-corner face index for multi-tile grids; entries are
    ``None`` for single-tile grids (including fold grids).

    ARGUMENTS
    ---------
    grid : `xgcm.Grid` instance
    mask : `xr.DataArray` instance of type `bool`

    RETURNS
    -------
    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list
    """
    arcs, closed, nf, Nyc, Nxc = _trace_and_drop(grid, mask)
    if get_facedim(grid) is not None:
        return _multitile_boundaries_from_mask(grid, arcs, closed, nf, Nyc, Nxc)
    return _single_tile_boundaries_from_mask(grid, arcs, closed)


def _remap_contour(c, o):
    """Map one `contourpy` polyline (center coordinates, first==last) to closed
    cell-corner index arrays `(i_c, j_c)`, with the corner-position offset `o`
    (1 for 'outer'/'left', 0 for 'right'; see `corner_offset`)."""
    i_c, j_c = c[:-1, 0], c[:-1, 1]
    i_n, j_n = i_c.copy(), j_c.copy()
    i_inc = np.roll(i_c, -1) - i_c
    j_inc = np.roll(j_c, -1) - j_c
    i_n[(i_c % 1) == 0.0] = (i_c - (i_inc < 0))[(i_c % 1) == 0.0] + o
    j_n[(j_c % 1) == 0.0] = (j_c - (j_inc < 0))[(j_c % 1) == 0.0] + o
    i_n[(i_c % 1) == 0.5] = np.floor(i_c[(i_c % 1) == 0.5]) + o
    j_n[(j_c % 1) == 0.5] = np.floor(j_c[(j_c % 1) == 0.5]) + o
    return loop(i_n).astype(np.int64), loop(j_n).astype(np.int64)


def _pad_center(grid, da):
    """Pad a center-point field one cell in each direction using the grid's own
    `xgcm` topology, so that only genuine topological seams carry real neighbours.

    Periodic axes, the north fold, and multi-tile `face_connections` supply real
    across-seam halo cells; every other boundary (walls, ``'extend'``, ...) is padded
    with NaN. Coercing non-seam boundaries to ``'fill'`` matters: ``'extend'`` would
    otherwise *replicate* the edge cell, so a wall segment would see an in-mask
    "neighbour" and be wrongly dropped as an interior seam. Mirrors
    `sectionate.gridutils.build_neighbor_maps`."""
    def _seam_or_fill(b):
        return b if (b == "periodic" or _is_fold_boundary(b)) else "fill"
    boundary = {ax: _seam_or_fill(grid.axes[ax].boundary) for ax in grid.axes}
    boundary_width = {ax: (1, 1) for ax in grid.axes}
    return pad(da, grid, boundary_width, boundary=boundary, fill_value=np.nan)


# Cells separated by a directed corner segment (ig,jg)->(ig+di,jg+dj), in the
# padded-center index frame where corner (jg,ig) straddles cells [jg:jg+2, ig:ig+2].
_SEG_CELLS = {
    (1, 0):  ((1, 1), (0, 1)),   # +i: north, south
    (-1, 0): ((1, 0), (0, 0)),   # -i: north, south
    (0, 1):  ((1, 1), (1, 0)),   # +j: east, west
    (0, -1): ((0, 1), (0, 0)),   # -j: east, west
}


def _trace_and_drop(grid, mask):
    """Shared front-end for both back-ends.

    `contourpy`-trace the mask per face (a single synthetic face for single-tile
    grids) and drop every boundary segment that lies on a seam between two in-mask
    cells -- i.e. interior to the region -- using the topology-aware halo `_pad_center`.
    Returns ``(arcs, closed, nf, Nyc, Nxc)`` in the padded-center corner-index frame,
    with each corner tagged ``(f, jg, ig)``:

    - ``arcs``   : open runs of kept segments (their endpoints sit on a seam);
    - ``closed`` : fully-kept contours as open-cyclic loops (N distinct corners,
      ``first != last``) -- back-ends finalize them.
    """
    facedim = get_facedim(grid)
    cdict = coord_dict(grid)
    Xc, Yc = cdict["X"]["center"], cdict["Y"]["center"]
    o = 1 - corner_offset(grid)

    if facedim is None:
        mask = mask.transpose(Yc, Xc)
        Nyc, Nxc = mask.shape
        nf = 1
        m = mask.values[None]
        Mpad = _pad_center(grid, mask.astype(float)).transpose(Yc, Xc).values[None]
    else:
        mask = mask.transpose(facedim, Yc, Xc)
        nf, Nyc, Nxc = mask.shape
        m = mask.values
        Mpad = _pad_center(grid, mask.astype(float)).transpose(facedim, ..., Yc, Xc).values

    arcs, closed = [], []
    for f in range(nf):
        z = np.pad(m[f].astype(float), 1)
        cs = contourpy.contour_generator(
            np.arange(-1, Nxc + 1), np.arange(-1, Nyc + 1), z
        ).create_contour(0.5)
        for c in cs:
            ig, jg = _remap_contour(c, o)
            ig, jg = ig[:-1], jg[:-1]  # open cyclic sequence
            N = len(ig)
            keep = np.ones(N, bool)
            for k in range(N):
                k2 = (k + 1) % N
                (aj, ai), (bj, bi) = _SEG_CELLS[(int(ig[k2] - ig[k]), int(jg[k2] - jg[k]))]
                # a segment is INTERNAL (a seam face between two in-mask cells) iff both
                # cells it separates are in-mask in the topology-aware halo
                if (Mpad[f, jg[k] + aj, ig[k] + ai] == 1.0
                        and Mpad[f, jg[k] + bj, ig[k] + bi] == 1.0):
                    keep[k] = False
            if keep.all():
                closed.append([(f, int(jg[k]), int(ig[k])) for k in range(N)])
                continue
            cut = np.where(~keep)[0]
            start = (cut[-1] + 1) % N
            run = []
            for t in range(N):
                k = (start + t) % N
                if keep[k]:
                    if not run:
                        run = [(f, int(jg[k]), int(ig[k]))]
                    run.append((f, int(jg[(k + 1) % N]), int(ig[(k + 1) % N])))
                elif run:
                    arcs.append(run)
                    run = []
            if run:
                arcs.append(run)
    return arcs, closed, nf, Nyc, Nxc


def _single_tile_boundaries_from_mask(grid, arcs, closed):
    """Single-tile back-end: stitch arcs into closed loops by physical coincidence at
    the seam. This handles walls, periodic axes, and the bipolar north fold uniformly:
    a region touching no seam traces exactly as a plain `contourpy` contour, while a
    region wrapping the periodic-X seam or straddling the fold is stitched into one
    seam-consistent loop whose crossings are expressed through *coincident* seam corners
    -- the form `sectionate` collapses to zero-length (dropped) faces. `f_c` is
    ``None`` (single tile)."""
    cdict = coord_dict(grid)
    Xq, Yq = cdict["X"]["corner"], cdict["Y"]["corner"]
    geo = get_geo_corners(grid)
    lon_c = geo["X"].transpose(Yq, Xq).values
    lat_c = geo["Y"].transpose(Yq, Xq).values

    def lockey(node):
        # physical position on the unit sphere: robust to longitude wrap (a coincident
        # seam corner may read 180 vs -180) and to the pole. Matches the physical
        # coincidence `sectionate` itself uses to collapse zero-length seam faces.
        _, jq, iq = node
        la, lo = np.deg2rad(float(lat_c[jq, iq])), np.deg2rad(float(lon_c[jq, iq]))
        return (round(np.cos(la) * np.cos(lo), 9),
                round(np.cos(la) * np.sin(lo), 9),
                round(np.sin(la), 9))

    ends = {}
    for ai, arc in enumerate(arcs):
        ends.setdefault(lockey(arc[0]), []).append((ai, True))
        ends.setdefault(lockey(arc[-1]), []).append((ai, False))
    used = [False] * len(arcs)
    loops = [lp + [lp[0]] for lp in closed]  # close the open-cyclic no-seam contours
    for a0 in range(len(arcs)):
        if used[a0]:
            continue
        lp, ai, at_start = [], a0, True
        while not used[ai]:
            used[ai] = True
            seg = arcs[ai] if at_start else arcs[ai][::-1]
            lp.extend(seg)  # keep BOTH coincident seam corners at each junction
            nxt = [(a, w) for (a, w) in ends.get(lockey(seg[-1]), []) if not used[a]]
            if not nxt:
                break
            ai, at_start = nxt[0]
        if lp[0] != lp[-1]:
            lp.append(lp[0])  # close (a coincident/zero-length edge sectionate drops)
        loops.append(lp)

    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list = [], [], [], [], []
    for lp in loops:
        j_c = np.array([n[1] for n in lp], dtype=np.int64)
        i_c = np.array([n[2] for n in lp], dtype=np.int64)
        i_c_list.append(i_c)
        j_c_list.append(j_c)
        f_c_list.append(None)
        lons_c_list.append(np.array([float(lon_c[n[1], n[2]]) for n in lp[:-1]]))
        lats_c_list.append(np.array([float(lat_c[n[1], n[2]]) for n in lp[:-1]]))
    return i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list


def _multitile_boundaries_from_mask(grid, arcs, closed, nf, Nyc, Nxc):
    """Multi-tile back-end: stitch arcs into loops across tile seams using the grid's
    `face_connections` topology, returning a per-corner face index `f_c`. Handles
    rotated/reversed seams whose two sides do not share coordinates (so coincidence
    stitching would fail) and produces the `f_c` sectionate's multi-tile transport
    needs."""
    facedim = get_facedim(grid)
    cdict = coord_dict(grid)
    Xc, Yc = cdict["X"]["center"], cdict["Y"]["center"]
    Xq, Yq = cdict["X"]["corner"], cdict["Y"]["corner"]
    geo = get_geo_corners(grid)
    lon_c = geo["X"].transpose(facedim, Yq, Xq).values
    lat_c = geo["Y"].transpose(facedim, Yq, Xq).values
    maps = build_neighbor_maps(grid, geo)

    cid = xr.DataArray(
        np.arange(nf * Nyc * Nxc, dtype=float).reshape(nf, Nyc, Nxc), dims=(facedim, Yc, Xc)
    )
    Cpad = _pad_center(grid, cid).transpose(facedim, ..., Yc, Xc).values

    def cellset(f, jg, ig):
        v = (Cpad[f, jg, ig], Cpad[f, jg, ig + 1], Cpad[f, jg + 1, ig], Cpad[f, jg + 1, ig + 1])
        return frozenset(None if np.isnan(x) else int(x) for x in v)

    def neighbours(f, j, i):
        out = []
        for d in NEIGHBOR_DIRECTIONS:
            fm, jm, im = maps[d]
            out.append((int(fm[f, j, i]), int(jm[f, j, i]), int(im[f, j, i])))
        return out

    # --- Stitch arcs into face-local loops by cell-set at endpoints ---
    ends = {}
    for ai, arc in enumerate(arcs):
        for node in (arc[0], arc[-1]):
            ends.setdefault(cellset(*node), []).append(ai)
    used = [False] * len(arcs)
    facelocal = list(closed)
    for a0 in range(len(arcs)):
        if used[a0]:
            continue
        lp, ai, from_start = [], a0, True
        while not used[ai]:
            used[ai] = True
            seg = arcs[ai] if from_start else arcs[ai][::-1]
            lp.extend(seg[:-1])
            tail = seg[-1]
            ks = cellset(*tail)
            nxt = [a for a in ends.get(ks, []) if not used[a]]
            if not nxt:
                lp.append(tail)
                break
            ai = nxt[0]
            from_start = cellset(*arcs[ai][0]) == ks
        facelocal.append(lp)

    # --- Convert face-local corners to native (f, j, i), grid-adjacent ---
    # native corner-array shape ('outer' has Nc+1 corners, 'left'/'right' have Nc)
    Nyq, Nxq = lon_c.shape[1], lon_c.shape[2]
    seed = {}
    for f in range(nf):
        for jn in range(Nyq):
            for inx in range(Nxq):
                seed.setdefault(cellset(f, jn, inx), (f, jn, inx))

    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list = [], [], [], [], []
    for lp in facelocal:
        targets = [cellset(*c) for c in lp]
        prev = seed.get(targets[0])
        if prev is None:
            continue
        nat = [prev]
        for k in range(1, len(targets)):
            cands = [prev] + neighbours(*prev)
            match = [c for c in cands if cellset(c[0], c[1], c[2]) == targets[k]]
            prev = match[0] if match else seed.get(targets[k], prev)
            nat.append(prev)

        # Repair seam crossings the cell-set match over-merged: where consecutive
        # corners are not grid-adjacent, insert the corner where they meet (the
        # neighbour of A that lies on B's tile and neighbours B).
        rep = []
        for k in range(len(nat)):
            a = nat[k]
            rep.append(a)
            b = nat[(k + 1) % len(nat)]
            if a != b and b not in neighbours(*a):
                bridge = [c for c in neighbours(*a) if c[0] == b[0] and c in neighbours(*b)]
                if bridge:
                    rep.append(bridge[0])

        seq = rep if rep[-1] == rep[0] else rep + [rep[0]]   # close exactly once
        f_c = np.array([c[0] for c in seq], dtype=np.int64)
        j_c = np.array([c[1] for c in seq], dtype=np.int64)
        i_c = np.array([c[2] for c in seq], dtype=np.int64)
        i_c_list.append(i_c)
        j_c_list.append(j_c)
        f_c_list.append(f_c)
        lons_c_list.append(np.array([float(lon_c[c[0], c[1], c[2]]) for c in seq[:-1]]))
        lats_c_list.append(np.array([float(lat_c[c[0], c[1], c[2]]) for c in seq[:-1]]))

    return i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list
