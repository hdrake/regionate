"""Optional end-to-end checks against real ocean-model grids.

These are skipped by default. Enable them by setting REGIONATE_REALDATA_TESTS=1
with the corresponding example file present locally under ``data/`` -- the MOM6
global example, and (for the ECCO checks) the LLC90 geometry downloaded by
``examples/load_example_ECCO_grid.py``. They validate the multi-model goal of the
overhaul on genuine MOM6 (tripolar) and ECCO (lat-lon-cap) output. The default
suite never reaches out to the network.
"""

import os
import numpy as np
import xarray as xr
import pytest

REALDATA = os.environ.get("REGIONATE_REALDATA_TESTS") == "1"
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data")
MOM6_FILE = os.path.join(
    DATA_DIR, "MOM6_global_example_vertically_integrated_mass_budget_v0_0_6.nc"
)
MOM6_HEAT_FILE = os.path.join(
    DATA_DIR, "MOM6_global_example_vertically_integrated_heat_budget_v0_0_6.nc"
)
ECCO_FILE = os.path.join(DATA_DIR, "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc")
EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "examples")


@pytest.mark.skipif(
    not (REALDATA and os.path.isfile(MOM6_FILE)),
    reason="set REGIONATE_REALDATA_TESTS=1 and provide the MOM6 global example file",
)
def test_mom6_global_box_mask_boundary_consistency():
    """On the real global MOM6 grid, a mid-latitude box's mask splits into connected
    components whose OWN traced boundaries close the discrete divergence theorem: a
    synthetic flux field's convergence over each component's cells equals the net flux
    through that component's boundary loops. This is the substantive consistency check
    (the boundaries really do enclose those cells) -- not the tautology that the
    component masks tile the input mask, which the labeling guarantees by construction."""
    import xgcm
    import sectionate as sec
    from regionate import GriddedRegion, MaskRegions

    ds = xr.open_dataset(MOM6_FILE).fillna(0.)
    grid = xgcm.Grid(
        ds,
        coords={"X": {"center": "xh", "outer": "xq"},
                "Y": {"center": "yh", "outer": "yq"}},
        padding={"X": "periodic", "Y": "extend"},
        autoparse_metadata=False,
    )
    lons = np.array([-40., -10., -10., -40.])
    lats = np.array([10., 10., 40., 40.])
    region = GriddedRegion("box", lons, lats, grid)
    assert int(region.mask.sum()) > 0

    mr = MaskRegions(region.mask, grid)
    assert len(mr.region_dict) >= 1
    # cheap sanity: the component masks tile the input mask exactly
    union = None
    for r in mr.region_dict.values():
        union = r.mask if union is None else (union | r.mask)
    assert bool((union.values == region.mask.values).all())

    # substance: convergence over each component's cells == flux through its boundaries
    ny, nx = ds.sizes["yh"], ds.sizes["xh"]
    umo = xr.DataArray(np.sin(0.011 * np.add.outer(np.arange(ny), np.arange(nx + 1))) + 0.3,
                       dims=("yh", "xq"))
    vmo = xr.DataArray(np.cos(0.013 * np.add.outer(np.arange(ny + 1), np.arange(nx))) - 0.2,
                       dims=("yq", "xh"))
    ds["umo"], ds["vmo"] = umo, vmo
    conv = -(grid.diff(ds["umo"], "X") + grid.diff(ds["vmo"], "Y"))
    for r in mr.region_dict.values():
        interior = float(conv.where(r.mask, 0.).sum())
        flux = 0.0
        for b in r.boundaries:
            t = sec.convergent_transport(grid, b.i_c, b.j_c, f_c=b.f_c,
                                         utr="umo", vtr="vmo", layer=None,
                                         positive_in=r.mask)
            flux += float(t["conv_mass_transport"].sum())
        assert np.isclose(flux, interior, rtol=1e-9, atol=1e-6)


@pytest.mark.skipif(
    not (REALDATA and os.path.isfile(MOM6_HEAT_FILE)),
    reason="set REGIONATE_REALDATA_TESTS=1 and provide the MOM6 heat-budget example file",
)
def test_mom6_arctic_fold_region_is_single_and_closes_budget():
    """On the real global MOM6 tripolar grid, the high-Arctic cold region straddles
    the bipolar north fold. With the fold boundary it must trace into a SINGLE region
    (rather than the two halves the ``Y='extend'`` wall gives), that region's boundary
    must reach the fold seam, and its advective-heat convergence must close the budget
    against the volume-integrated tendency (the discrete divergence theorem). This is
    the worked example of ``examples/3_Arctic_heat_CM4p25.ipynb`` as a regression test,
    and exercises the real-grid seam coincidence (antimeridian/pole) the fold stitch
    relies on."""
    import xgcm
    import sectionate as sec
    from regionate import MaskRegions

    ds = xr.open_dataset(MOM6_HEAT_FILE).fillna(0.)
    ds = ds.expand_dims(["z_l"]).assign_coords({
        "z_l": xr.DataArray([3000], dims=("z_l",)),
        "z_i": xr.DataArray([0, 6000], dims=("z_i",))})
    coords = {"X": {"center": "xh", "outer": "xq"},
              "Y": {"center": "yh", "outer": "yq"}}
    mask = ((ds["tos"].squeeze() < 0.) & (ds["geolat"] > 0)).compute()
    mask = xr.DataArray(mask.values, dims=("yh", "xh"),
                        coords={"geolon": ds.geolon, "geolat": ds.geolat})
    seam = ds.sizes["yq"] - 1

    def nbnd(r):
        return sum(len(np.asarray(b.i_c)) for b in r.boundaries)

    def touches_seam(r):
        return any((np.asarray(b.j_c) == seam).any() for b in r.boundaries)

    def large(grid):
        return [r for r in MaskRegions(mask, grid).region_dict.values() if nbnd(r) > 100]

    grid_extend = xgcm.Grid(ds, coords=coords, padding={"X": "periodic", "Y": "extend"},
                            autoparse_metadata=False)
    grid_fold = xgcm.Grid(ds, coords=coords, padding={"X": "periodic", "Y": {"fold": "corner"}},
                          autoparse_metadata=False)

    big_extend = large(grid_extend)
    big_fold = large(grid_fold)
    # the fold merges the two seam-split halves: fewer large seam-touching regions
    seam_extend = sum(touches_seam(r) for r in big_extend)
    seam_fold = sum(touches_seam(r) for r in big_fold)
    assert seam_fold < seam_extend
    # the single biggest fold region spans the seam (traces across the fold)
    biggest = max(big_fold, key=nbnd)
    assert touches_seam(biggest)

    # discrete divergence theorem: boundary heat convergence == volume tendency
    regions = MaskRegions(mask, grid_fold).region_dict
    total = 0.0
    for r in regions.values():
        for b in r.boundaries:
            dsec = sec.convergent_transport(
                grid_fold, b.i_c, b.j_c, f_c=b.f_c, utr="T_adx", vtr="T_ady",
                layer="z_l", interface="z_i", outname="cht", positive_in=r.mask)
            total += float(dsec["cht"].sum("z_l").isel(time=0).sum(["sect"]).values)
    dheatdt = (ds["T_advection_xy"] * ds["areacello"]).sum("z_l")
    tend = float(dheatdt.where(mask).sum(["xh", "yh"]).isel(time=0).values)
    assert np.isclose(total, tend, rtol=1e-3)


# Opt-in gate for the ECCO end-to-end tests, mirroring the MOM6 test above: they
# require REGIONATE_REALDATA_TESTS=1 *and* a locally-present LLC90 geometry file, so
# the default suite (and CI) never downloads from the network. NB this must decorate
# the test functions themselves -- a skipif marker on the `_load_ecco` helper is inert
# (pytest only honours markers on collected `test_*` functions).
_requires_ecco = pytest.mark.skipif(
    not (REALDATA and os.path.isfile(ECCO_FILE)),
    reason="set REGIONATE_REALDATA_TESTS=1 and download the ECCO LLC90 geometry "
           "(see examples/load_example_ECCO_grid.py)",
)


def _load_ecco():
    import sys
    sys.path.insert(0, os.path.abspath(EXAMPLES_DIR))
    from load_example_ECCO_grid import load_ECCO_LLC90_grid, atlantic_basin_mask
    grid = load_ECCO_LLC90_grid(data_dir=os.path.abspath(DATA_DIR))
    return grid, atlantic_basin_mask


@_requires_ecco
def test_ecco_llc90_seam_region_stitches_across_tiles():
    """On the real ECCOv4r4 lat-lon-cap (LLC90, native MITgcm 'left') grid, a
    contiguous region straddling the tile-1/tile-2 seam is traced as a single
    boundary loop whose per-corner face index spans both tiles."""
    from regionate import MaskRegions
    grid, _ = _load_ecco()
    lon, lat = grid._ds["geolon"], grid._ds["geolat"]
    mask = ((lon > -30) & (lon < 20) & (lat > -5) & (lat < 25)).compute()

    regions = MaskRegions(mask, grid).region_dict
    assert len(regions) == 1
    faces = set().union(*(set(np.asarray(b.f_c).tolist())
                          for b in regions[0].boundaries))
    assert faces == {1, 2}


@_requires_ecco
def test_ecco_atlantic_basin_boundary_and_transports():
    """The published Atlantic basin (regionmask) on the LLC90 grid traces as a
    boundary spanning many tiles across rotated seams and the Arctic cap, and --
    critically -- the boundary is grid-adjacent everywhere, so it converts to
    velocity (u,v) faces via sectionate (the flux-divergence consistency that is
    the whole point of the package)."""
    import sectionate as sec
    from regionate import MaskRegions
    grid, atlantic_basin_mask = _load_ecco()

    mask = atlantic_basin_mask(grid)
    regions = MaskRegions(mask, grid).region_dict
    basin = max(regions.values(), key=lambda r: int(r.mask.sum()))
    faces = set().union(*(set(np.asarray(b.f_c).tolist()) for b in basin.boundaries))
    assert len(faces) >= 4                      # spans many tiles (rotated seams)

    # Every boundary loop must be convertible to velocity faces -- this is the
    # strong test: it requires every consecutive corner pair to be grid-adjacent.
    for b in basin.boundaries:
        lons_uv, lats_uv = sec.uvcoords_from_qindices(
            grid, b.i_c, b.j_c, f_c=b.f_c)
        assert len(lons_uv) > 0


def _ecco_synthetic_uv(grid):
    """An arbitrary synthetic transport field on the native LLC90 staggering.

    The southern boundary fold (the j_g=0 rows of tiles 0 and 3 fold back onto
    themselves and each other) stores each physical edge twice, so a transport
    field is only self-consistent there if the two storages are antisymmetric.
    MITgcm guarantees this trivially -- those edges are domain walls carrying no
    flux -- so the synthetic field zeroes them the same way."""
    nf, ny, nx = (grid._ds.sizes[d] for d in ("tile", "j", "i"))
    g = np.arange(nf * ny * nx, dtype=float).reshape(nf, ny, nx)
    umo = xr.DataArray(np.sin(g * 0.013) + 0.3, dims=("tile", "j", "i_g"))
    vmo = xr.DataArray(np.cos(g * 0.017) - 0.2, dims=("tile", "j_g", "i"))
    vmo[{"tile": [0, 3], "j_g": 0}] = 0.0
    return umo, vmo


def _ecco_convergence(grid, umo, vmo):
    """Cell convergence from the native transports via the outer corner topology.

    ``sectionate.gridutils.outer_topology(grid).padded_transports`` resolves each
    face's missing edge slots to the *stored* velocity of that physical edge, so
    the convergence telescopes exactly: its global sum is identically zero and it
    is exactly consistent with boundary fluxes read from the same native arrays.
    (``padded_transports`` is used rather than xgcm's ``grid.diff(...,
    other_component=...)`` because it is xgcm-independent and resolves edges
    stored on no face -- walls, the lon=-115 cut, the 4th Arctic vertex -- which
    no halo pad can supply a value for. xgcm#749 fixes the bare-``DataArray``
    pad path across rotated/reversed seams; the dict form was always exact.)"""
    from sectionate.gridutils import outer_topology
    Uo, Vo = outer_topology(grid).padded_transports(
        umo.transpose("tile", "j", "i_g"), vmo.transpose("tile", "j_g", "i")
    )
    conv = (Uo[:, :, :-1] - Uo[:, :, 1:]) + (Vo[:, :-1, :] - Vo[:, 1:, :])
    return xr.DataArray(conv, dims=("tile", "j", "i"))


def _boundary_flux(grid, mask, umo, vmo):
    """Net flux into `mask` through its traced boundary loops, via sectionate."""
    import sectionate as sec
    from regionate import MaskRegions

    ds = grid._ds
    total = 0.0
    ds["umo"], ds["vmo"] = umo, vmo
    for r in MaskRegions(mask, grid).region_dict.values():
        for b in r.boundaries:
            t = sec.convergent_transport(
                grid, b.i_c, b.j_c, b.f_c, utr="umo", vtr="vmo",
                layer=None, positive_in=r.mask,
            )
            total += float(t["conv_mass_transport"].sum())
    return total


@_requires_ecco
def test_ecco_atlantic_basin_obeys_discrete_divergence_theorem():
    """The whole point of the package: a region's budget must close against the
    fluxes through its traced boundary. On the real LLC90 grid, for the full
    Atlantic basin (spanning rotated seams), the net flux through every boundary
    velocity face (summed over all loops) equals the flux convergence summed over
    the masked cells -- exactly, for an arbitrary transport field."""
    grid, atlantic_basin_mask = _load_ecco()
    mask = atlantic_basin_mask(grid)
    umo, vmo = _ecco_synthetic_uv(grid)
    conv = _ecco_convergence(grid, umo, vmo)
    interior = float(conv.where(mask, 0.).sum())
    flux = _boundary_flux(grid, mask, umo, vmo)
    assert np.isclose(flux, interior, rtol=1e-12, atol=1e-6)


@_requires_ecco
@pytest.mark.parametrize("case", ["south_cap", "north_cap", "vertex_annulus_rot",
                                  "vertex_annulus_latlon", "latlon_rot_box"])
def test_ecco_hard_topology_regions_close_exactly(case):
    """Exact mask<->boundary closure on the LLC90 grid's hardest topology: the
    south-pole boundary fold and lon=-115 grid cut (south cap), the Arctic cap
    with its rotated seams (north cap), nested annuli around a rotated-rotated
    and a lat-lon<->rotated 4-face cube-vertex junction, and a box crossing a
    lat-lon<->rotated seam."""
    grid, _ = _load_ecco()
    lon, lat = grid._ds["geolon"], grid._ds["geolat"]

    def geodist(lon0, lat0):
        la0, lo0 = np.deg2rad(lat0), np.deg2rad(lon0)
        la, lo = np.deg2rad(lat), np.deg2rad(lon)
        return np.rad2deg(np.arccos(np.clip(
            np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0),
            -1., 1.)))

    masks = {
        "south_cap": lat < -60.,
        "north_cap": lat > 70.,
        "vertex_annulus_rot": (geodist(-128., 9.97) < 3.5) & (geodist(-128., 9.97) > 1.4),
        "vertex_annulus_latlon": (geodist(-38., 9.97) < 3.5) & (geodist(-38., 9.97) > 1.4),
        "latlon_rot_box": (((lon - (-60.)) % 360.) <= 50.) & (lat > 20.) & (lat < 50.),
    }
    mask = masks[case].compute()
    umo, vmo = _ecco_synthetic_uv(grid)
    conv = _ecco_convergence(grid, umo, vmo)
    interior = float(conv.where(mask, 0.).sum())
    flux = _boundary_flux(grid, mask, umo, vmo)
    assert np.isclose(flux, interior, rtol=1e-12, atol=1e-6)
