"""Boundary tracing and stitching across tile seams on multi-tile grids
(grids defined by xgcm `face_connections`, e.g. lat-lon-cap / cubed-sphere)."""

import numpy as np
import xarray as xr
import xgcm
import pytest

import sectionate as sec
from regionate.boundaries import grid_boundaries_from_mask
from regionate import MaskRegions


def rotated_two_tile_grid(Nc=4):
    """A native 'left'-staggered 2-tile grid whose tiles meet at a ROTATED seam:
    face 0's +X edge connects to face 1's Y axis (an X->Y, 90-degree connection).
    Each tile is given its own disjoint coordinate block (face 1 offset by +100),
    so the seam's two corner representations never coincide in lon/lat -- a region
    spanning the seam can therefore only be stitched via the grid topology
    (cell adjacency), not coordinate coincidence."""
    ng = Nc  # native 'left' corners are the same size as centers
    LONc = np.zeros((2, ng, ng)); LATc = np.zeros((2, ng, ng))
    LON = np.zeros((2, Nc, Nc));  LAT = np.zeros((2, Nc, Nc))
    for f in range(2):
        off = 100 * f
        LONc[f] = off + np.arange(ng)[None, :]
        LATc[f] = off + np.arange(ng)[:, None]
        LON[f] = off + np.arange(Nc)[None, :] + 0.5
        LAT[f] = off + np.arange(Nc)[:, None] + 0.5
    ds = xr.Dataset(coords={
        "i": ("i", np.arange(Nc)), "j": ("j", np.arange(Nc)),
        "i_g": ("i_g", np.arange(ng)), "j_g": ("j_g", np.arange(ng)),
        "face": ("face", [0, 1]),
        "geolon": (("face", "j", "i"), LON), "geolat": (("face", "j", "i"), LAT),
        "geolon_c": (("face", "j_g", "i_g"), LONc),
        "geolat_c": (("face", "j_g", "i_g"), LATc),
    })
    fc = {"face": {0: {"X": (None, (1, "Y", False))},
                   1: {"Y": ((0, "X", False), None)}}}
    grid_left = xgcm.Grid(
        ds, coords={"X": {"center": "i", "left": "i_g"},
                    "Y": {"center": "j", "left": "j_g"}},
        padding="fill", fill_value=np.nan,
        face_connections=fc, autoparse_metadata=False,
    )
    return grid_left


def test_rotated_seam_region_stitches_into_one_loop():
    grid = rotated_two_tile_grid(Nc=4)
    # A region spanning the rotated seam: the seam between the tiles is internal
    # and the region's boundary must be a single loop spanning both faces --
    # stitched across the rotated seam by topology alone (coords don't coincide).
    # The mask stays clear of the tiles' open outer walls, whose corner points a
    # native 'left' grid does not store (see `test_open_wall_boundary_raises`).
    arr = np.zeros(grid._ds["geolon"].shape, dtype=bool)
    # clear of every corner the tiling does not store: face 0's north wall row
    # and its j=0 row (whose SE corner is the seam's wall end), face 1's east
    # and north wall rows.
    arr[0, 1:3, 2:4] = True
    arr[1, 0:3, 0:3] = True
    mask = xr.DataArray(arr, dims=grid._ds["geolon"].dims,
                        coords=grid._ds["geolon"].coords)
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1
    assert set(np.asarray(f_l[0]).tolist()) == {0, 1}


def test_open_wall_boundary_raises():
    """A mask reaching an open wall whose corner points are stored on no face
    (the high-side rows/columns of a native 'left' tile with no neighbour there)
    cannot be expressed in native (i_c, j_c, f_c) indices; regionate must say so
    rather than emit fabricated corners."""
    grid = rotated_two_tile_grid(Nc=4)
    mask = xr.ones_like(grid._ds["geolon"]).astype(bool)
    with pytest.raises((ValueError, RuntimeError), match="Mask boundary"):
        grid_boundaries_from_mask(grid, mask)


def two_face_grid(Nc=3):
    """Two faces side-by-side in longitude: face 0 spans [0, 90], face 1 [90, 180],
    joined at a single seam (face 0's right X edge -> face 1's left X edge). Carries
    both tracer-center (geolon/geolat) and corner (geolon_c/geolat_c) coordinates."""
    ng = Nc + 1
    yq = np.linspace(-45, 45, ng)
    yh = 0.5 * (yq[:-1] + yq[1:])
    lonq = [np.linspace(0, 90, ng), np.linspace(90, 180, ng)]
    lonh = [0.5 * (l[:-1] + l[1:]) for l in lonq]

    LONc = np.stack([np.broadcast_to(lonq[f], (ng, ng)) for f in range(2)])
    LATc = np.stack([np.broadcast_to(yq[:, None], (ng, ng)) for f in range(2)])
    LON = np.stack([np.broadcast_to(lonh[f], (Nc, Nc)) for f in range(2)])
    LAT = np.stack([np.broadcast_to(yh[:, None], (Nc, Nc)) for f in range(2)])

    ds = xr.Dataset(
        {},
        coords={
            "xq": (("xq",), np.arange(ng)), "yq": (("yq",), np.arange(ng)),
            "xh": (("xh",), np.arange(Nc)), "yh": (("yh",), np.arange(Nc)),
            "face": (("face",), [0, 1]),
            "geolon_c": (("face", "yq", "xq"), LONc),
            "geolat_c": (("face", "yq", "xq"), LATc),
            "geolon": (("face", "yh", "xh"), LON),
            "geolat": (("face", "yh", "xh"), LAT),
        },
    )
    fc = {"face": {0: {"X": (None, (1, "X", False))},
                   1: {"X": ((0, "X", False), None)}}}
    return xgcm.Grid(
        ds,
        coords={"X": {"outer": "xq", "center": "xh"},
                "Y": {"outer": "yq", "center": "yh"}},
        padding="fill", fill_value=np.nan,
        face_connections=fc, autoparse_metadata=False,
    )


def make_mask(grid, cells):
    """cells: dict mapping face index -> list of (j, i) center cells to set True."""
    arr = np.zeros_like(grid._ds.geolon.values, dtype=bool)
    for f, lst in cells.items():
        for (j, i) in lst:
            arr[f, j, i] = True
    return xr.DataArray(arr, dims=grid._ds.geolon.dims, coords=grid._ds.geolon.coords)


def test_boundary_obeys_discrete_divergence_theorem():
    """A traced multi-tile boundary must obey the discrete divergence theorem:
    for any flux field, the net flux through the boundary's velocity faces equals
    the flux convergence summed over the masked cells. This is the property that
    makes regionate budgets consistent. Uses a region spanning the (non-rotated)
    tile seam and a seam-consistent synthetic transport field.

    (This synthetic grid is symmetric 'outer', so the shared seam U-face is stored
    on BOTH tiles and must be made single-valued by hand; on a native 'left' grid --
    e.g. real ECCO -- the seam face is stored once and any transport field works, as
    ``test_ecco_atlantic_basin_obeys_discrete_divergence_theorem`` checks.)"""
    grid = two_face_grid(Nc=6)
    Nc = grid._ds.sizes["xh"]; ng = Nc + 1
    # synthetic face transports; the shared seam U-face (face0 xq=Nc == face1 xq=0)
    # must be single-valued for the flux field to be physically consistent.
    umo = np.sin(np.arange(2 * Nc * ng).reshape(2, Nc, ng) * 0.07) + 0.3
    umo[1, :, 0] = umo[0, :, Nc]
    umo = xr.DataArray(umo, dims=("face", "yh", "xq"))
    vmo = xr.DataArray(np.cos(np.arange(2 * ng * Nc).reshape(2, ng, Nc) * 0.05) - 0.2,
                       dims=("face", "yq", "xh"))

    mask = make_mask(grid, {0: [(j, i) for j in range(1, 5) for i in range(3, Nc)],
                            1: [(j, i) for j in range(1, 5) for i in range(0, 3)]})
    convergence = float((-(grid.diff(umo, "X") + grid.diff(vmo, "Y"))).where(mask, 0.).sum())

    U, V = umo.values, vmo.values
    flux = 0.0
    for region in MaskRegions(mask, grid).region_dict.values():
        for loop_ in region.boundaries:
            uv = sec.uvindices_from_qindices(grid, loop_.i_c, loop_.j_c, f_c=loop_.f_c)
            for k in range(len(uv["var"])):
                if uv["var"][k] == "0":
                    continue
                f, i, j = int(uv["face"][k]), int(uv["i"][k]), int(uv["j"][k])
                flux += int(uv["Lsign"][k]) * (U[f, j, i] if uv["var"][k] == "U" else V[f, j, i])

    assert np.isclose(convergence, flux, atol=1e-9)


def test_seam_spanning_region_stitches_into_one_loop():
    grid = two_face_grid(Nc=3)
    # East column of face 0 + west column of face 1: a strip straddling the seam.
    mask = make_mask(grid, {0: [(0, 2), (1, 2), (2, 2)],
                            1: [(0, 0), (1, 0), (2, 0)]})
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)

    assert len(i_l) == 1                                  # a single stitched loop
    faces = set(np.asarray(f_l[0]).tolist())
    assert faces == {0, 1}                                # spans both faces
    # Every physical corner appears exactly once (a shared seam corner is not
    # duplicated), so every boundary edge is a real velocity face: the 3x2-cell
    # strip has 10 boundary corners and 10 velocity faces.
    lon_uv, lat_uv = sec.uvcoords_from_qindices(grid, i_l[0], j_l[0], f_c=f_l[0])
    assert len(lon_uv) == 10
    assert len(lon_l[0]) == 10


def test_seam_terminating_region_keeps_seam_edge():
    grid = two_face_grid(Nc=3)
    # East column of face 0 only: its east neighbour (face 1) is outside the mask,
    # so the seam edge at lon=90 is a real boundary and must be kept.
    mask = make_mask(grid, {0: [(0, 2), (1, 2), (2, 2)]})
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)

    assert len(i_l) == 1
    assert set(np.asarray(f_l[0]).tolist()) == {0}       # never leaves face 0
    assert np.isclose(np.asarray(lon_l[0]), 90.).any()   # seam edge retained
    # No internal seam faces here: every corner contributes a velocity face.
    lon_uv, _ = sec.uvcoords_from_qindices(grid, i_l[0], j_l[0], f_c=f_l[0])
    assert len(lon_uv) == len(lon_l[0])


def test_interior_cell_single_face_box():
    grid = two_face_grid(Nc=3)
    mask = make_mask(grid, {1: [(1, 1)]})
    i_l, j_l, f_l, lon_l, lat_l = grid_boundaries_from_mask(grid, mask)
    assert len(i_l) == 1
    assert set(np.asarray(f_l[0]).tolist()) == {1}
    assert len(lon_l[0]) == 4                             # a 4-corner cell box


def test_maskregions_threads_face_index():
    grid = two_face_grid(Nc=3)
    mask = make_mask(grid, {0: [(0, 2), (1, 2), (2, 2)],
                            1: [(0, 0), (1, 0), (2, 0)]})
    regions = MaskRegions(mask, grid).region_dict
    assert len(regions) == 1
    region = regions[0]
    assert len(region.boundaries) == 1
    loop_ = region.boundaries[0]
    assert loop_.f_c is not None
    assert set(np.asarray(loop_.f_c).tolist()) == {0, 1}


def test_boundary_to_mask_rasterizes_across_tiles():
    """boundary -> mask on a multi-tile grid: rasterizing a lon/lat polygon that
    straddles the tile seam must fill cells on BOTH faces. regionmask only accepts
    1D/2D lon/lat, so this exercises `mask_from_grid_boundaries`' per-tile
    rasterize-and-stitch (`rasterize_per_tile`) rather than handing it the 3D
    (face, y, x) coordinate array."""
    from regionate.grid_conform import mask_from_grid_boundaries
    grid = two_face_grid(Nc=6)  # face 0 spans lon [0,90], face 1 [90,180]
    # A box straddling the lon=90 seam, interior to the tiles' outer walls.
    lons = np.array([45., 135., 135., 45.])
    lats = np.array([-20., -20., 20., 20.])
    mask = mask_from_grid_boundaries(lons, lats, grid)

    assert "face" in mask.dims                              # face dimension retained
    assert bool(mask.sel(face=0).any())                    # filled on face 0 ...
    assert bool(mask.sel(face=1).any())                    # ... and face 1
    # Exactly the cells whose centers fall inside the lon/lat box (centers are
    # strictly interior to the box edges, so the closed polygon selects them).
    lon = grid._ds.geolon.values
    lat = grid._ds.geolat.values
    expected = (lon > 45.) & (lon < 135.) & (lat > -20.) & (lat < 20.)
    assert np.array_equal(mask.values, expected)


def _geodist(lon, lat, lon0, lat0):
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    la0, lo0 = np.deg2rad(lat0), np.deg2rad(lon0)
    return np.rad2deg(np.arccos(np.clip(
        np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0), -1., 1.)))


@pytest.mark.parametrize("lon0,lat0,radius,faces_expected", [
    (0., 45., 18., {0, 4}),       # one rotated seam
    (45., 40., 16., {0, 1, 4}),   # three faces meeting near a cube corner
    (90., 45., 18., {1, 4}),      # another rotated seam
])
def test_cube_rotated_seam_obeys_divergence_theorem(lon0, lat0, radius, faces_expected):
    """A region straddling a ROTATED tile seam must obey the discrete divergence
    theorem, checked against an INDEPENDENT ground truth: xgcm's own convergence
    (`grid.diff` with `other_component`) equals the net flux through regionate's
    traced boundary. This closes the rotated-seam coverage gap that otherwise runs
    only in the gated real-ECCO tests.

    The fixture is a physically-valid cubed-sphere (`cube_grid.cube_left_grid`),
    whose seam corner coordinates coincide and whose rotated `face_connections`
    make `xgcm.diff` a legitimate oracle -- unlike `rotated_two_tile_grid`, whose
    deliberately non-physical offset coords make it a stitching-only fixture on
    which no divergence test (xgcm.diff or padded_transports) closes. These regions
    stay clear of the cube vertices that live on no face -- NOT to dodge them, but
    because `xgcm.diff` is not a valid oracle there (its halo fabricates a value at
    an edge stored on no face). Those hardest points are covered by
    `test_cube_vertex_junction_closes_for_nondivergent_flow`, whose div-free oracle
    IS valid at the vertices."""
    from cube_grid import cube_left_grid, Nc
    grid, _ = cube_left_grid()
    lon = grid._ds["geolon"].values
    lat = grid._ds["geolat"].values
    m = _geodist(lon, lat, lon0, lat0) < radius
    mask = xr.DataArray(m, dims=grid._ds["geolon"].dims, coords=grid._ds["geolon"].coords)
    assert set(np.where(m.any(axis=(1, 2)))[0].tolist()) == faces_expected  # straddles the seam(s)

    rng = np.random.default_rng(1)
    umo = xr.DataArray(rng.standard_normal((6, Nc, Nc)), dims=("face", "j", "i_g"))
    vmo = xr.DataArray(rng.standard_normal((6, Nc, Nc)), dims=("face", "j_g", "i"))

    # independent ground truth: xgcm's own convergence across the rotated seam
    conv = -(grid.diff(umo, "X", other_component={"Y": vmo})
             + grid.diff(vmo, "Y", other_component={"X": umo}))
    interior = float(conv.where(mask, 0.).sum())

    # net flux through the traced boundary, via sectionate
    U, V = umo.values, vmo.values
    i_l, j_l, f_l, _, _ = grid_boundaries_from_mask(grid, mask)
    flux = 0.0
    for k in range(len(i_l)):
        uv = sec.uvindices_from_qindices(grid, i_l[k], j_l[k], f_c=f_l[k])
        for t in range(len(uv["var"])):
            if uv["var"][t] == "0":
                continue
            f, i, j = int(uv["face"][t]), int(uv["i"][t]), int(uv["j"][t])
            flux += int(uv["Lsign"][t]) * (U[f, j, i] if uv["var"][t] == "U" else V[f, j, i])

    assert np.isclose(interior, flux, atol=1e-9)


def test_cube_vertex_junction_closes_for_nondivergent_flow():
    """The hardest topology: a region whose boundary wraps a cube-vertex junction --
    three faces and their rotated seams meeting at a point stored on at most one
    face, INCLUDING the two vertices stored on NO face (the cube analogue of LLC90's
    Arctic vertex). `xgcm.diff` cannot be the oracle here (its halo fabricates a
    value at an edge stored on no face), so closure is checked against the
    topology-independent fact that a closed loop's net transport of a non-divergent
    (streamfunction) flow is exactly zero -- which holds even at the vertices."""
    from cube_grid import cube_left_grid
    from sectionate.topology import corner_topology as outer_topology
    grid, _ = cube_left_grid()
    ot = outer_topology(grid)
    deg = np.array([len(a) for a in ot.node_adj])
    native = ot.node_native[:, 0] >= 0
    assert np.count_nonzero(~native) == 2   # the fixture really has 2 unstored vertices
    # both unstored vertices (on no face) + one stored 3-valent vertex
    vertices = list(np.where(~native)[0]) + [int(np.where((deg == 3) & native)[0][0])]

    lon = grid._ds["geolon"].values
    lat = grid._ds["geolat"].values
    u, v = grid._ds["u"].values, grid._ds["v"].values   # non-divergent streamfunction flow
    for n in vertices:
        lo0, la0 = float(ot.node_lon[n]), float(ot.node_lat[n])
        m = _geodist(lon, lat, lo0, la0) < 25.
        mask = xr.DataArray(m, dims=grid._ds["geolon"].dims, coords=grid._ds["geolon"].coords)
        assert len(set(np.where(m.any(axis=(1, 2)))[0].tolist())) >= 3   # truly wraps the junction

        i_l, j_l, f_l, _, _ = grid_boundaries_from_mask(grid, mask)
        flux = 0.0
        for k in range(len(i_l)):
            uv = sec.uvindices_from_qindices(grid, i_l[k], j_l[k], f_c=f_l[k])
            for t in range(len(uv["var"])):
                if uv["var"][t] == "0":
                    continue
                f, i, j = int(uv["face"][t]), int(uv["i"][t]), int(uv["j"][t])
                flux += int(uv["Lsign"][t]) * (u[f, j, i] if uv["var"][t] == "U" else v[f, j, i])
        assert np.isclose(flux, 0., atol=1e-9)
