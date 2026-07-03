import contourpy
import numpy as np
import xarray as xr
from xgcm.padding import pad, _is_fold_boundary

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

    The tracer follows the grid's own topology, so it is correct for every
    `xgcm.Grid` supported by the fold/`face_connections` padding:

    - Plain single-tile grids (periodic and/or walled) are traced with `contourpy`.
    - Single-tile grids carrying a bipolar/tripolar north **fold**
      (`boundary={..., "Y": {"fold": ...}}`) are traced with `contourpy`, but
      segments lying on the fold seam between two in-mask fold-neighbour cells are
      dropped (they are interior, not boundary) and the surviving arcs are stitched
      into loops by physical coincidence at the seam -- so a fold-straddling region
      yields a single loop with no spurious seam boundary.
    - Multi-tile grids (`face_connections`, e.g. lat-lon-cap or cubed-sphere) are
      traced face-by-face and stitched across tile seams using the grid topology.

    In every case a region yields closed boundary loops whose corners map cleanly to
    velocity faces via `sectionate.uvindices_from_qindices`, so integrating a flux
    over the boundary reproduces its convergence over the masked cells exactly.

    Returns lists with a common length equal to the number of discrete boundary
    loops. `f_c_list` holds the per-corner face index for multi-tile grids; its
    entries are `None` for single-tile grids (including fold grids).

    ARGUMENTS
    ---------
    grid : `xgcm.Grid` instance
    mask : `xr.DataArray` instance of type `bool`

    RETURNS
    -------
    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list
    """
    if get_facedim(grid) is not None:
        return _multitile_boundaries_from_mask(grid, mask)
    if _has_fold(grid):
        return _fold_boundaries_from_mask(grid, mask)
    return _single_tile_boundaries_from_mask(grid, mask)


def _has_fold(grid):
    """True if any axis carries a bipolar/tripolar north-fold boundary (xgcm#711)."""
    for ax in grid.axes:
        try:
            if _is_fold_boundary(grid.axes[ax].boundary):
                return True
        except (AttributeError, TypeError):
            pass
    return False


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


def _single_tile_boundaries_from_mask(grid, mask):
    """Trace mask boundaries on a single-tile grid with `contourpy`."""
    o = 1 - corner_offset(grid)
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
        i_c_new, j_c_new = _remap_contour(c, o)
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


def _fold_boundaries_from_mask(grid, mask):
    """Trace mask boundaries on a single-tile grid carrying a bipolar north fold
    (xgcm#711 ``boundary={..., "Y": {"fold": ...}}``).

    Like the plain single-tile path this uses ``contourpy``, but it additionally
    (1) drops boundary segments that lie on the fold seam between two in-mask
    fold-neighbour cells -- they are interior to the region, not boundary -- using
    the grid's own fold-aware halo (`_pad_center`), and (2) stitches the surviving
    arcs into closed loops by physical coincidence at the seam. A region straddling
    the fold therefore yields a single loop whose fold crossings are expressed
    through *coincident* seam corners -- the representation `sectionate` collapses to
    zero-length (dropped) faces -- so the boundary's velocity faces reproduce the
    region's flux convergence exactly. `f_c` is ``None`` (single tile).
    """
    cdict = coord_dict(grid)
    Xc, Yc = cdict["X"]["center"], cdict["Y"]["center"]
    Xq, Yq = cdict["X"]["corner"], cdict["Y"]["corner"]
    geo = get_geo_corners(grid)
    lon_c = geo["X"].transpose(Yq, Xq).values
    lat_c = geo["Y"].transpose(Yq, Xq).values
    o = 1 - corner_offset(grid)

    mask = mask.transpose(Yc, Xc)
    Nyc, Nxc = mask.shape
    m = mask.values
    # fold-aware halo: a top-edge cell's north neighbour is its across-fold mirror
    Mpad = _pad_center(grid, mask.astype(float)).transpose(Yc, Xc).values

    def lockey(jq, iq):
        # physical position on the unit sphere: robust to longitude wrap (a coincident
        # seam corner may read 180 vs -180) and to the pole. Matches the physical
        # coincidence `sectionate` itself uses to collapse zero-length seam faces.
        la, lo = np.deg2rad(float(lat_c[jq, iq])), np.deg2rad(float(lon_c[jq, iq]))
        return (round(np.cos(la) * np.cos(lo), 9),
                round(np.cos(la) * np.sin(lo), 9),
                round(np.sin(la), 9))

    # --- Stage 1: contourpy -> arcs of KEEP segments (+ already-closed loops) ---
    z = np.pad(m.astype(float), 1)
    contours = contourpy.contour_generator(
        np.arange(-1, Nxc + 1), np.arange(-1, Nyc + 1), z
    ).create_contour(0.5)
    arcs, closed = [], []
    for c in contours:
        ig, jg = _remap_contour(c, o)
        ig, jg = ig[:-1], jg[:-1]  # open cyclic sequence
        N = len(ig)
        keep = np.ones(N, bool)
        for k in range(N):
            k2 = (k + 1) % N
            (aj, ai), (bj, bi) = _SEG_CELLS[(int(ig[k2] - ig[k]), int(jg[k2] - jg[k]))]
            # a segment is INTERNAL (a fold-seam face between two in-mask cells) iff
            # both cells it separates are in-mask in the fold-aware halo
            if (Mpad[jg[k] + aj, ig[k] + ai] == 1.0
                    and Mpad[jg[k] + bj, ig[k] + bi] == 1.0):
                keep[k] = False
        if keep.all():
            closed.append([(int(jg[k]), int(ig[k])) for k in range(N)]
                          + [(int(jg[0]), int(ig[0]))])
            continue
        cut = np.where(~keep)[0]
        start = (cut[-1] + 1) % N
        run = []
        for t in range(N):
            k = (start + t) % N
            if keep[k]:
                if not run:
                    run = [(int(jg[k]), int(ig[k]))]
                run.append((int(jg[(k + 1) % N]), int(ig[(k + 1) % N])))
            elif run:
                arcs.append(run)
                run = []
        if run:
            arcs.append(run)

    # --- Stage 2: stitch arcs into loops by endpoint coincidence (seam collapses) ---
    ends = {}
    for ai, arc in enumerate(arcs):
        ends.setdefault(lockey(*arc[0]), []).append((ai, True))
        ends.setdefault(lockey(*arc[-1]), []).append((ai, False))
    used = [False] * len(arcs)
    loops = list(closed)
    for a0 in range(len(arcs)):
        if used[a0]:
            continue
        lp, ai, at_start = [], a0, True
        while not used[ai]:
            used[ai] = True
            seg = arcs[ai] if at_start else arcs[ai][::-1]
            lp.extend(seg)  # keep BOTH coincident seam corners at each junction
            nxt = [(a, w) for (a, w) in ends.get(lockey(*seg[-1]), []) if not used[a]]
            if not nxt:
                break
            ai, at_start = nxt[0]
        if lp[0] != lp[-1]:
            lp.append(lp[0])  # close (a coincident/zero-length edge sectionate drops)
        loops.append(lp)

    # --- Stage 3: emit index arrays (closed) + coordinates (open) ---
    i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list = [], [], [], [], []
    for lp in loops:
        j_c = np.array([c[0] for c in lp], dtype=np.int64)
        i_c = np.array([c[1] for c in lp], dtype=np.int64)
        i_c_list.append(i_c)
        j_c_list.append(j_c)
        f_c_list.append(None)
        lons_c_list.append(np.array([float(lon_c[c[0], c[1]]) for c in lp[:-1]]))
        lats_c_list.append(np.array([float(lat_c[c[0], c[1]]) for c in lp[:-1]]))
    return i_c_list, j_c_list, f_c_list, lons_c_list, lats_c_list


def _pad_center(grid, da):
    """Pad a center-point field one cell in each direction using the grid's own
    `xgcm` topology (periodic / `face_connections` / fold). Genuine walls are
    filled with NaN. Mirrors `sectionate.gridutils.build_neighbor_maps`."""
    boundary = {ax: grid.axes[ax].boundary for ax in grid.axes}
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


def _multitile_boundaries_from_mask(grid, mask):
    """Trace and stitch mask boundaries across tile seams on a multi-tile grid."""
    facedim = get_facedim(grid)
    cdict = coord_dict(grid)
    Xc, Yc = cdict["X"]["center"], cdict["Y"]["center"]
    Xq, Yq = cdict["X"]["corner"], cdict["Y"]["corner"]
    geo = get_geo_corners(grid)
    lon_c = geo["X"].transpose(facedim, Yq, Xq).values
    lat_c = geo["Y"].transpose(facedim, Yq, Xq).values
    maps = build_neighbor_maps(grid, geo)

    mask = mask.transpose(facedim, Yc, Xc)
    nf, Nyc, Nxc = mask.shape
    Mpad = _pad_center(grid, mask.astype(float)).transpose(facedim, ..., Yc, Xc).values
    cid = xr.DataArray(
        np.arange(nf * Nyc * Nxc, dtype=float).reshape(nf, Nyc, Nxc), dims=(facedim, Yc, Xc)
    )
    Cpad = _pad_center(grid, cid).transpose(facedim, ..., Yc, Xc).values
    m = mask.values

    def cellset(f, jg, ig):
        v = (Cpad[f, jg, ig], Cpad[f, jg, ig + 1], Cpad[f, jg + 1, ig], Cpad[f, jg + 1, ig + 1])
        return frozenset(None if np.isnan(x) else int(x) for x in v)

    def neighbours(f, j, i):
        out = []
        for d in NEIGHBOR_DIRECTIONS:
            fm, jm, im = maps[d]
            out.append((int(fm[f, j, i]), int(jm[f, j, i]), int(im[f, j, i])))
        return out

    # --- Stage 1+2: contourpy per face -> open arcs of KEEP segments + closed loops ---
    arcs, closed = [], []
    for f in range(nf):
        z = np.pad(m[f].astype(float), 1)
        cs = contourpy.contour_generator(
            np.arange(-1, Nxc + 1), np.arange(-1, Nyc + 1), z
        ).create_contour(0.5)
        for c in cs:
            ig, jg = _remap_contour(c, 1)
            ig, jg = ig[:-1], jg[:-1]  # open cyclic sequence
            N = len(ig)
            keep = np.ones(N, bool)
            for k in range(N):
                k2 = (k + 1) % N
                (aj, ai), (bj, bi) = _SEG_CELLS[(int(ig[k2] - ig[k]), int(jg[k2] - jg[k]))]
                # a segment is INTERNAL (cut) iff both cells it separates are in-mask
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

    # --- Stage 3: stitch arcs into face-local loops by cell-set at endpoints ---
    ends = {}
    for ai, arc in enumerate(arcs):
        for f, jg, ig in (arc[0], arc[-1]):
            ends.setdefault(cellset(f, jg, ig), []).append(ai)
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
