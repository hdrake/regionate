import contourpy
import numpy as np
import xarray as xr
# scipy is a hard dependency of `sectionate` (this package's core sibling), so it is
# always available; used for C-speed seam-aware connected-component labeling.
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components as _connected_components_csr
from xgcm.padding import pad
try:
    # north-fold boundary detector; only present in xgcm with north-fold support (xgcm#711)
    from xgcm.padding import _is_fold_padding as _is_fold_boundary
except ImportError:  # pragma: no cover - fall back so `import regionate` works on any xgcm
    from collections.abc import Mapping

    def _is_fold_boundary(boundary):
        return isinstance(boundary, Mapping) and "fold" in boundary

from .utilities import loop
from sectionate.gridutils import (
    get_facedim,
    coord_dict,
)
from sectionate.topology import corner_topology


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
    arcs are then stitched into closed loops on the grid's corner topology
    (`sectionate.topology.corner_topology`): every traced corner resolves to a physical
    corner *node*, arcs join by node identity -- uniform across a periodic wrap, a
    bipolar fold, rotated and reversed tile seams, and cube-vertex junctions -- and each
    node is emitted as the native corner that stores it, with its face index `f_c`.

    One back-end, for every grid. A single-tile grid used to be stitched separately, by
    matching rounded positions on the unit sphere, which meant this package and
    `sectionate` had to agree on a coincidence tolerance and did not; it also merged
    corners a grid distinguishes but happens to place at one point, which is what a
    tripolar cap's singular meridian does to a whole column of them.

    In every case the returned loops' corners map cleanly to velocity faces via
    `sectionate.uvindices_from_qindices`, so integrating a flux over the boundary
    reproduces its convergence over the masked cells exactly (a property of the full
    loop *set* -- an annulus, fold, or multi-tile region may return several loops).

    Returns lists with a common length equal to the number of discrete boundary loops.
    `f_c_list` holds the per-corner face index; a single-tile grid is one face, so it
    is zeros there rather than ``None``.

    ARGUMENTS
    ---------
    grid : `xgcm.Grid` instance
    mask : `xr.DataArray` instance of type `bool`

    RETURNS
    -------
    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list
    """
    arcs, closed, nf, Nyc, Nxc = _trace_and_drop(grid, mask)
    return _boundaries_from_arcs(grid, arcs, closed, nf, Nyc, Nxc)


def _remap_contour(c, o):
    """Map one `contourpy` polyline (center coordinates, first==last) to closed
    cell-corner index arrays `(i_c, j_c)`, with the corner offset `o`
    (always 1 here: the corner topology is indexed on the 'outer' lattice)."""
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
    "neighbour" and be wrongly dropped as an interior seam.

    **This diverges from `sectionate.topology.CornerTopology` and that matters.**
    It reads the grid's own metadata only, so it does not see identifications
    declared through `sectionate.topology.declare_identifications` -- LLC90's
    southern boundary fold, for one. Cells either side of such a seam look
    unconnected here, so `connected_components` splits a component that spans it.
    The boundary tracer is unaffected, because it annihilates the doubly-traced
    edges within a single trace, so budgets still close; it is region *identity*
    that is wrong. On LLC90 the cells concerned are Antarctic land, which is why
    nothing observes it yet. The fix is to take cell adjacency from the corner
    topology, which has seen those identifications."""
    def _seam_or_fill(b):
        return b if (b == "periodic" or _is_fold_boundary(b)) else "fill"
    # Only pad axes whose dimensions `da` actually carries (issue #24): a center-point
    # field spans just the horizontal (X/Y) tracer dims, but a 3D `grid` also declares a
    # Z axis. Padding a Z axis absent from `da.dims` makes `xgcm.pad` raise
    # `KeyError: None of the DataArray's dims (...) were found in axis coords`.
    axes = [ax for ax in grid.axes
            if set(grid.axes[ax].coords.values()) & set(da.dims)]
    padding = {ax: _seam_or_fill(grid.axes[ax].padding) for ax in axes}
    padding_width = {ax: (1, 1) for ax in axes}
    return pad(da, grid, padding_width, padding=padding, fill_value=np.nan)


def connected_components(grid, mask):
    """Label a boolean cell `mask` into connected components on the grid's own
    seam-aware 4-adjacency.

    Two in-mask cells belong to the same component iff they are edge-neighbours
    *through the grid's topology*: periodic axes, the north fold, and multi-tile
    `face_connections` all supply real across-seam neighbours, exactly as they do
    for the boundary tracer -- both read the same topology-aware halo,
    `_pad_center`. A planar labeller (e.g. ``scipy.ndimage.label``) would instead
    wrongly split a component that wraps a seam and merge cells that are only
    planar-adjacent across a wall, so we label on the padded halo instead.

    The id field padded here holds each cell's *global index*, so every padded
    halo slot carries the real across-seam neighbour's id (and NaN at a genuine
    wall). Two in-mask cells adjacent through the halo are then unioned; component
    indices are assigned in first-appearance (flat) order and are otherwise not
    significant.

    ARGUMENTS
    ---------
    grid : `xgcm.Grid` instance
    mask : `xr.DataArray` of bool over the grid's tracer-center dims

    RETURNS
    -------
    labels : `xr.DataArray` of int over the same dims as `mask` (transposed to the
        grid's center dims); each in-mask cell holds its component index in
        ``range(ncomp)`` and every out-of-mask cell holds ``-1``.
    ncomp : int
        The number of connected components (0 for an all-False mask).
    """
    facedim = get_facedim(grid)
    cdict = coord_dict(grid)
    Xc, Yc = cdict["X"]["center"], cdict["Y"]["center"]

    if facedim is None:
        m = mask.transpose(Yc, Xc).astype(bool)
        mvals = m.values[None]  # synthetic single face
    else:
        m = mask.transpose(facedim, Yc, Xc).astype(bool)
        mvals = m.values
    nf, Nyc, Nxc = mvals.shape
    size = nf * Nyc * Nxc

    # A center-point id field: each cell holds its own global index. Padding it with
    # the grid's topology puts the *neighbour's* id in each halo slot (NaN at walls).
    ids = np.arange(size, dtype=float).reshape(nf, Nyc, Nxc)
    ids_da = xr.DataArray(ids.reshape(m.shape), dims=m.dims, coords=m.coords)
    if facedim is None:
        Ip = _pad_center(grid, ids_da).transpose(Yc, Xc).values[None]
    else:
        Ip = _pad_center(grid, ids_da).transpose(facedim, Yc, Xc).values

    mflat = mvals.reshape(-1)
    self_id = np.arange(size, dtype=np.int64)
    # native cell (f,j,i) sits at Ip[f, j+1, i+1]; its four edge-neighbours:
    neighbours = (
        Ip[:, 1:-1, 2:],   # east
        Ip[:, 1:-1, :-2],  # west
        Ip[:, 2:, 1:-1],   # north
        Ip[:, :-2, 1:-1],  # south
    )

    # Collect the seam-aware edges between in-mask cells, then label with scipy's
    # C-speed connected-components. scipy is a hard dependency of `sectionate` (this
    # package's core sibling), so this adds no new dependency. We build the graph on
    # only the in-mask cells (remapped to a compact index) so out-of-mask cells do
    # not each become their own component.
    in_idx = np.flatnonzero(mflat)              # global ids of in-mask cells (sorted)
    n_in = in_idx.size
    pos = np.full(size, -1, dtype=np.int64)     # global id -> compact in-mask index
    pos[in_idx] = np.arange(n_in)

    src_parts, dst_parts = [], []
    for nbr in neighbours:
        n = nbr.reshape(-1)
        valid = mflat & ~np.isnan(n)            # self in-mask, neighbour exists
        s_v = self_id[valid]
        n_v = n[valid].astype(np.int64)
        both = mflat[n_v]                        # neighbour also in-mask
        src_parts.append(pos[s_v[both]])
        dst_parts.append(pos[n_v[both]])
    src = np.concatenate(src_parts) if src_parts else np.empty(0, np.int64)
    dst = np.concatenate(dst_parts) if dst_parts else np.empty(0, np.int64)

    graph = coo_matrix(
        (np.ones(src.size, dtype=np.int8), (src, dst)), shape=(n_in, n_in)
    ).tocsr()
    ncomp, comp = _connected_components_csr(graph, directed=False)

    # Relabel so component ids appear in first-in-mask-cell (flat) order -- stable and
    # independent of scipy's internal labeling.
    _, first = np.unique(comp, return_index=True)
    remap = np.empty(ncomp, dtype=np.int64)
    remap[comp[np.sort(first)]] = np.arange(ncomp)
    comp = remap[comp]

    labels_flat = np.full(size, -1, dtype=np.int64)
    labels_flat[in_idx] = comp

    labels_arr = labels_flat.reshape(nf, Nyc, Nxc)
    if facedim is None:
        labels_arr = labels_arr[0]
    labels = xr.DataArray(labels_arr, dims=m.dims, coords=m.coords)
    return labels, ncomp


# Cells separated by a directed corner segment (ig,jg)->(ig+di,jg+dj), in the
# padded-center index frame where corner (jg,ig) straddles cells [jg:jg+2, ig:ig+2].
_SEG_CELLS = {
    (1, 0):  ((1, 1), (0, 1)),   # +i: north, south
    (-1, 0): ((1, 0), (0, 0)),   # -i: north, south
    (0, 1):  ((1, 1), (1, 0)),   # +j: east, west
    (0, -1): ((0, 1), (0, 0)),   # -j: east, west
}


def _trace_and_drop(grid, mask):
    """Front-end for `_boundaries_from_arcs`.

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
    # Corner indices in the 'outer'-lattice frame (`o = 1`, corner k between cells
    # k-1 and k), which is what the corner topology's node grid is indexed by --
    # on every grid, whatever its native staggering.
    o = 1

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



def _boundaries_from_arcs(grid, arcs, closed, nf, Nyc, Nxc):
    """Stitch the traced arcs into closed loops on the grid's corner topology.

    One back-end for every grid. Every traced corner -- given by `_trace_and_drop`
    in the 'outer'-lattice frame ``(f, jg, ig)`` -- resolves to a physical corner
    *node*, and that is uniform across a periodic wrap, a bipolar fold, rotated or
    reversed tile seams, and cube-vertex junctions alike. There is nothing left for
    a single-tile grid to do differently: it used to stitch by matching rounded
    positions on the unit sphere, which meant regionate and sectionate had to agree
    on a coincidence tolerance, and which quietly merged corners a grid distinguishes
    but places at one point -- a tripolar cap's singular meridian. Stitching is now:

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
    ot = corner_topology(grid)
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
                f"any face (corner slot {tuple(int(x) for x in ot.reps_of(lp[k])[0])}"
                + (f", near lon={ot.node_lon[lp[k]]:.2f}, lat={ot.node_lat[lp[k]]:.2f}"
                   if ot.node_position_known[lp[k]]
                   else ", whose position the grid does not determine either")
                + "); it cannot be expressed in native "
                "(i_c, j_c, f_c) indices."
            )
        f_c_list.append(nat[:, 0].astype(np.int64))
        j_c_list.append(nat[:, 1].astype(np.int64))
        i_c_list.append(nat[:, 2].astype(np.int64))
        lons_c_list.append(ot.node_lon[lp[:-1]].astype(float))
        lats_c_list.append(ot.node_lat[lp[:-1]].astype(float))

    return i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list
