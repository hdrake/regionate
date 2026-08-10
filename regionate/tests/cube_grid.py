"""A synthetic cubed-sphere in native 'left' (MITgcm) staggering — a *valid*
rotated multi-tile fixture for closure tests.

Unlike the deliberately non-physical `rotated_two_tile_grid` (offset coords, only
usable for stitching tests), this cube is generated geometry-first (gnomonic
faces), its `face_connections` are *derived* from corner coincidence, and its
corners carry real physical lon/lat that coincide across every seam. Its rotated
seams are therefore physically consistent, so `xgcm.diff` is a legitimate
independent oracle on them (away from the two cube vertices that live on no face).

Adapted from `sectionate/tests/test_cube_left_grid.py` (only the grid builder is
copied here; the streamfunction flow is available for div-free checks).
"""

import numpy as np
import xarray as xr
import xgcm

Nc = 8

# cube face frames: (center, local +x unit, local +y unit) of each face.
_BASE_FRAMES = [
    ((1, 0, 0), (0, 1, 0), (0, 0, 1)),    # +X
    ((0, 1, 0), (-1, 0, 0), (0, 0, 1)),   # +Y
    ((-1, 0, 0), (0, -1, 0), (0, 0, 1)),  # -X
    ((0, -1, 0), (1, 0, 0), (0, 0, 1)),   # -Y
    ((0, 0, 1), (0, 1, 0), (-1, 0, 0)),   # +Z (top)
    ((0, 0, -1), (0, 1, 0), (1, 0, 0)),   # -Z (bottom)
]


def _rotated_frames(rots):
    frames = []
    for (c, u, v), k in zip(_BASE_FRAMES, rots):
        u, v = np.array(u, dtype=float), np.array(v, dtype=float)
        for _ in range(k):        # in-plane quarter turn (preserves handedness)
            u, v = v, -u
        frames.append((np.array(c, dtype=float), u, v))
    return frames


def _edge_endpoints(frames, f, which):
    c, u, v = frames[f]
    pts = {
        "Xlo": ((c - u - v), (c - u + v)),
        "Xhi": ((c + u - v), (c + u + v)),
        "Ylo": ((c - u - v), (c + u - v)),
        "Yhi": ((c - u + v), (c + u + v)),
    }[which]
    return tuple(tuple(np.round(p / np.linalg.norm(p), 9)) for p in pts)


def _shared_edge_sides(frames, f, f2):
    for w in ("Xlo", "Xhi", "Ylo", "Yhi"):
        e = set(_edge_endpoints(frames, f, w))
        for w2 in ("Xlo", "Xhi", "Ylo", "Yhi"):
            if e == set(_edge_endpoints(frames, f2, w2)):
                return (w[-2:], w2[-2:])
    return None


def _find_all_low_high_rotations():
    """Per-face quarter-turns making every cube edge a low<->high gluing (the
    orientation xgcm's `face_connections` padding handles reliably)."""
    import itertools
    for rots in itertools.product(range(4), repeat=6):
        frames = _rotated_frames(rots)
        ok = True
        for f in range(6):
            for f2 in range(f + 1, 6):
                sides = _shared_edge_sides(frames, f, f2)
                if sides is not None and sides[0] == sides[1]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return rots
    raise AssertionError("no all-low<->high cube orientation found")


_FRAMES = _rotated_frames(_find_all_low_high_rotations())


def _face_xyz(face, a, b):
    c, u, v = _FRAMES[face]
    p = c[:, None, None] + a[None] * u[:, None, None] + b[None] * v[:, None, None]
    return p / np.linalg.norm(p, axis=0)


def _outer_corners_xyz():
    e = np.linspace(-1.0, 1.0, Nc + 1)
    a, b = np.meshgrid(e, e, indexing="xy")
    out = np.empty((6, Nc + 1, Nc + 1, 3))
    for f in range(6):
        out[f] = np.moveaxis(_face_xyz(f, a, b), 0, -1)
    return out


def _lonlat(xyz):
    lon = np.rad2deg(np.arctan2(xyz[..., 1], xyz[..., 0]))
    lat = np.rad2deg(np.arcsin(np.clip(xyz[..., 2], -1.0, 1.0)))
    return lon, lat


def _derive_face_connections(cxyz):
    """Read the cube topology off the geometry, in xgcm's
    (neighbor_face, neighbor_axis, reverse) convention."""
    def edge_line(f, which):
        C = cxyz[f]
        return {"Xlo": C[:, 0], "Xhi": C[:, Nc], "Ylo": C[0, :], "Yhi": C[Nc, :]}[which]

    def key(p):
        return tuple(np.round(p, 9))

    lines = {(f, w): edge_line(f, w) for f in range(6)
             for w in ("Xlo", "Xhi", "Ylo", "Yhi")}
    conns = {f: {"X": [None, None], "Y": [None, None]} for f in range(6)}
    for (f, w), L in lines.items():
        match = None
        for (f2, w2), L2 in lines.items():
            if f2 == f:
                continue
            ends, ends2 = {key(L[0]), key(L[Nc])}, {key(L2[0]), key(L2[Nc])}
            if ends == ends2:
                match = (f2, w2)
                break
        assert match is not None, f"unmatched cube edge {(f, w)}"
        f2, w2 = match
        axis2 = "X" if w2 in ("Xlo", "Xhi") else "Y"
        side = 0 if w.endswith("lo") else 1
        side2 = 0 if w2.endswith("lo") else 1
        conns[f]["X" if w.startswith("X") else "Y"][side] = (f2, axis2, side == side2)
    return {"face": {f: {ax: tuple(v) for ax, v in d.items()}
                     for f, d in conns.items()}}


def _psi(xyz):
    """A smooth streamfunction of position (single-valued at every shared corner)."""
    return (1.3 * xyz[..., 0] - 0.7 * xyz[..., 1]
            + 0.4 * xyz[..., 2] + 0.9 * xyz[..., 0] * xyz[..., 1] * xyz[..., 2])


def cube_left_grid():
    """Native 'left'-staggered cubed-sphere grid carrying an exactly non-divergent
    streamfunction flow in ``u``/``v``. Returns ``(grid, psi)`` where ``psi`` is the
    corner streamfunction ``(6, Nc+1, Nc+1)``."""
    cxyz = _outer_corners_xyz()
    lon_o, lat_o = _lonlat(cxyz)
    psi = _psi(cxyz)

    lonc, latc = lon_o[:, :Nc, :Nc], lat_o[:, :Nc, :Nc]
    e = np.linspace(-1.0, 1.0, Nc + 1)
    m = 0.5 * (e[:-1] + e[1:])
    am, bm = np.meshgrid(m, m, indexing="xy")
    cen = np.stack([np.moveaxis(_face_xyz(f, am, bm), 0, -1) for f in range(6)])
    lonh, lath = _lonlat(cen)

    u = psi[:, 1:, :Nc] - psi[:, :-1, :Nc]         # (6, Nc, Nc) at (j, i_g)
    v = -(psi[:, :Nc, 1:] - psi[:, :Nc, :-1])      # (6, Nc, Nc) at (j_g, i)

    ds = xr.Dataset(
        {"u": (("face", "j", "i_g"), u), "v": (("face", "j_g", "i"), v)},
        coords={
            "i": ("i", np.arange(Nc)), "j": ("j", np.arange(Nc)),
            "i_g": ("i_g", np.arange(Nc)), "j_g": ("j_g", np.arange(Nc)),
            "face": ("face", np.arange(6)),
            "geolon": (("face", "j", "i"), lonh), "geolat": (("face", "j", "i"), lath),
            "geolon_c": (("face", "j_g", "i_g"), lonc),
            "geolat_c": (("face", "j_g", "i_g"), latc),
        },
    )
    fc = _derive_face_connections(cxyz)
    grid = xgcm.Grid(
        ds, coords={"X": {"center": "i", "left": "i_g"},
                    "Y": {"center": "j", "left": "j_g"}},
        padding="fill", fill_value=np.nan,
        face_connections=fc, autoparse_metadata=False,
    )
    return grid, psi
