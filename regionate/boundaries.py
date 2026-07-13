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
    outer_topology,
)


def grid_boundaries_from_mask(grid, mask):
    """Find the cell-corner boundaries that enclose a boolean cell `mask`.

    The tracer follows the grid's own topology, so it stitches a region across every
    seam its `xgcm.Grid` declares -- periodic axes, the bipolar/tripolar north
    **fold** (`padding={..., "Y": {"fold": ...}}`), and multi-tile `face_connections`
    (lat-lon-cap / cubed-sphere). A region that wraps a seam therefore yields a single
    seam-consistent boundary rather than pieces split at the seam with spurious seam
    faces.

    All cases share one front-end (`_trace_and_drop`): `contourpy` traces the mask per
    face, and boundary segments lying on a seam between two in-mask cells are dropped
    as interior (using the grid's own topology-aware halo, `_pad_center`). The surviving
    arcs are then stitched into closed loops by one of two back-ends:

    - single-tile (periodic and/or walled, incl. the fold) -- stitched by physical
      coincidence at the seam; face index `f_c` is ``None``;
    - multi-tile (`face_connections`) -- stitched on the grid's outer (shared-corner)
      lattice from `sectionate.gridutils.outer_topology`: every traced corner resolves
      to a physical corner *node*, arcs join by node identity (uniform across rotated
      and reversed seams, cube-vertex junctions, and grid cuts/folds), and each node is
      emitted as the native corner that stores it, with its face index `f_c`.

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
    padding = {ax: _seam_or_fill(grid.axes[ax].padding) for ax in grid.axes}
    padding_width = {ax: (1, 1) for ax in grid.axes}
    return pad(da, grid, padding_width, padding=padding, fill_value=np.nan)


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
    # Single-tile: corner indices in the native frame. Multi-tile: always in the
    # outer-lattice frame (`o=1`, corner k between cells k-1 and k), which is what
    # `outer_topology`'s node grid is indexed by, whatever the native staggering.
    o = 1 if facedim is not None else 1 - corner_offset(grid)

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
    """Multi-tile back-end: stitch the traced arcs into loops on the grid's outer
    (shared-corner) corner topology (`sectionate.gridutils.outer_topology`).

    Every traced corner -- given by `_trace_and_drop` in the outer-lattice frame
    ``(f, jg, ig)`` -- resolves to a physical corner *node*, which is uniform
    across rotated/reversed seams, 4-face cube-vertex junctions, the pole, and
    grid cuts/folds. Stitching is then simply:

    1. decompose arcs into directed corner-to-corner segments of node pairs
       (dropping zero-length segments between coincident corners, e.g. a fold
       pleat tip);
    2. remove *both* copies of any edge traced twice -- an edge is traced once
       per adjacent in-mask cell, so a double appearance means in-mask cells on
       both sides of an undeclared seam (the two coincident sides of a grid
       cut or boundary fold, e.g. under Antarctica on the LLC grid): interior,
       not boundary;
    3. chain the surviving fragments end-to-end by node identity into closed
       loops;
    4. emit each node as the native corner that stores it, `(i_c, j_c, f_c)`.
    """
    ot = outer_topology(grid)
    node_id, node_native = ot.node_id, ot.node_native

    def node_of(f, jg, ig):
        n = int(node_id[f, jg, ig])
        if n < 0:
            raise ValueError(
                f"Mask boundary passes through corner slot (face={f}, j={jg}, i={ig}) "
                "that could not be resolved to a physical grid corner."
            )
        return n

    # --- 1. directed segments as node pairs (arcs and fully-closed contours) ---
    fragments = []  # each: list of node ids, len >= 2
    for arc in arcs:
        fragments.append([node_of(*c) for c in arc])
    for lp in closed:
        seq = [node_of(*c) for c in lp]
        fragments.append(seq + [seq[0]])  # close the open-cyclic contour

    segments = []
    for frag in fragments:
        for a, b in zip(frag[:-1], frag[1:]):
            if a != b:  # coincident corners (zero-length edge) carry no boundary
                segments.append((a, b))

    # --- 2. annihilate edges traced from both sides: interior to a cut/fold ---
    count = {}
    for a, b in segments:
        key = (min(a, b), max(a, b))
        count[key] = count.get(key, 0) + 1
    if any(c > 2 for c in count.values()):
        raise RuntimeError(
            "A boundary edge was traced more than twice; the mask topology is "
            "inconsistent with the grid's corner topology."
        )
    kept = [(a, b) for (a, b) in segments if count[(min(a, b), max(a, b))] == 1]

    # --- 3. chain directed segments into closed loops by node identity ---
    # Each surviving directed segment is used exactly once. Segments inherit
    # contourpy's orientation (the in-mask side is consistently to one side),
    # so following out-segments from each end node reproduces closed loops.
    out_by_node = {}
    for k, (a, b) in enumerate(kept):
        out_by_node.setdefault(a, []).append(k)
    used = [False] * len(kept)
    loops = []
    for k0 in range(len(kept)):
        if used[k0]:
            continue
        used[k0] = True
        a0, b = kept[k0]
        lp = [a0, b]
        while b != a0:
            nxt = [k for k in out_by_node.get(b, []) if not used[k]]
            if not nxt:
                raise RuntimeError(
                    "Mask boundary does not close on the grid's corner topology "
                    f"(dead end at corner node {b})."
                )
            k = nxt[0]
            used[k] = True
            b = kept[k][1]
            lp.append(b)
        loops.append(lp)

    # --- 4. native corners and coordinates per node ---
    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list = [], [], [], [], []
    for lp in loops:
        nat = node_native[lp]
        if (nat[:, 0] < 0).any():
            k = int(np.where(nat[:, 0] < 0)[0][0])
            raise ValueError(
                "Mask boundary passes through a grid corner that is not stored on "
                f"any face (near lon={ot.node_lon[lp[k]]:.2f}, "
                f"lat={ot.node_lat[lp[k]]:.2f}); it cannot be expressed in native "
                "(i_c, j_c, f_c) indices."
            )
        f_c_list.append(nat[:, 0].astype(np.int64))
        j_c_list.append(nat[:, 1].astype(np.int64))
        i_c_list.append(nat[:, 2].astype(np.int64))
        lons_c_list.append(ot.node_lon[lp[:-1]].astype(float))
        lats_c_list.append(ot.node_lat[lp[:-1]].astype(float))

    return i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list
