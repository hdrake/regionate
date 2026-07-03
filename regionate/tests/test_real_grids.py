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
    """On the real global MOM6 grid, a mid-latitude box defined by its boundary
    yields a mask whose own traced boundary re-encloses the same mask."""
    import xgcm
    from regionate import GriddedRegion, MaskRegions

    ds = xr.open_dataset(MOM6_FILE).fillna(0.)
    grid = xgcm.Grid(
        ds,
        coords={"X": {"center": "xh", "outer": "xq"},
                "Y": {"center": "yh", "outer": "yq"}},
        boundary={"X": "periodic", "Y": "extend"},
        autoparse_metadata=False,
    )
    lons = np.array([-40., -10., -10., -40.])
    lats = np.array([10., 10., 40., 40.])
    region = GriddedRegion("box", lons, lats, grid)
    assert int(region.mask.sum()) > 0

    # Re-tracing the mask must recover a region enclosing exactly the same cells.
    retraced = MaskRegions(region.mask, grid).region_dict
    assert len(retraced) >= 1
    union = None
    for r in retraced.values():
        union = r.mask if union is None else (union | r.mask)
    assert bool((union.values == region.mask.values).all())


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

    def large(grid):
        return [r for r in MaskRegions(mask, grid).region_dict.values()
                if len(r.i_c) > 100]

    grid_extend = xgcm.Grid(ds, coords=coords, boundary={"X": "periodic", "Y": "extend"},
                            autoparse_metadata=False)
    grid_fold = xgcm.Grid(ds, coords=coords, boundary={"X": "periodic", "Y": {"fold": "corner"}},
                          autoparse_metadata=False)

    big_extend = large(grid_extend)
    big_fold = large(grid_fold)
    # the fold merges the two seam-split halves: fewer large seam-touching regions
    seam_extend = sum((np.asarray(r.j_c) == seam).any() for r in big_extend)
    seam_fold = sum((np.asarray(r.j_c) == seam).any() for r in big_fold)
    assert seam_fold < seam_extend
    # the single biggest fold region spans the seam (traces across the fold)
    biggest = max(big_fold, key=lambda r: len(r.i_c))
    assert (np.asarray(biggest.j_c) == seam).any()

    # discrete divergence theorem: boundary heat convergence == volume tendency
    regions = MaskRegions(mask, grid_fold).region_dict
    total = 0.0
    for r in regions.values():
        dsec = sec.convergent_transport(
            grid_fold, r.i_c, r.j_c, f_c=r.f_c, utr="T_adx", vtr="T_ady",
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
    assert set(np.asarray(regions[0].f_c).tolist()) == {1, 2}


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
    basin = max(regions.values(), key=lambda r: len(r.lons_c))
    faces = set(np.asarray(basin.f_c).tolist())
    assert len(faces) >= 4                      # spans many tiles (rotated seams)

    # The basin boundary must be convertible to velocity faces -- this is the
    # strong test: it requires every consecutive corner pair to be grid-adjacent.
    lons_uv, lats_uv = sec.uvcoords_from_qindices(
        grid, basin.i_c, basin.j_c, f_c=basin.f_c)
    assert len(lons_uv) > 0


@_requires_ecco
def test_ecco_atlantic_basin_obeys_discrete_divergence_theorem():
    """The whole point of the package: a region's budget must close against the
    fluxes through its traced boundary. On the real LLC90 grid, for the full
    Atlantic basin (spanning rotated seams), the net flux through every boundary
    velocity face (summed over all loops) equals the flux convergence summed over
    the masked cells -- to machine precision, for an arbitrary transport field.

    The cell-centred convergence is taken with xgcm's vector-aware ``grid.diff``
    (``other_component=``): across a 90-degree LLC seam the U-component rotates into
    the neighbour's V-component, so differencing the components as independent scalars
    would be wrong there (this is expected xgcm behaviour, not a bug -- see xgcm's
    vector-padding API)."""
    import sectionate as sec
    from regionate import MaskRegions
    grid, atlantic_basin_mask = _load_ecco()
    mask = atlantic_basin_mask(grid)

    nf, ny, nx = (grid._ds.sizes[d] for d in ("tile", "j", "i"))
    g = np.arange(nf * ny * nx, dtype=float).reshape(nf, ny, nx)
    umo = xr.DataArray(np.sin(g * 0.013) + 0.3, dims=("tile", "j", "i_g"))
    vmo = xr.DataArray(np.cos(g * 0.017) - 0.2, dims=("tile", "j_g", "i"))

    divU = grid.diff({"X": umo}, "X", other_component={"Y": vmo},
                     to="center", boundary="fill", fill_value=np.nan)
    divV = grid.diff({"Y": vmo}, "Y", other_component={"X": umo},
                     to="center", boundary="fill", fill_value=np.nan)
    convergence = float((-(divU + divV)).where(mask, 0.).sum())

    U, V = umo.transpose("tile", ...).values, vmo.transpose("tile", ...).values
    flux = 0.0
    for r in MaskRegions(mask, grid).region_dict.values():
        uv = sec.uvindices_from_qindices(grid, r.i_c, r.j_c, f_c=r.f_c)
        for k in range(len(uv["var"])):
            if uv["var"][k] == "0":
                continue
            f, i, j = int(uv["face"][k]), int(uv["i"][k]), int(uv["j"][k])
            flux += int(uv["Lsign"][k]) * (U[f, j, i] if uv["var"][k] == "U" else V[f, j, i])

    assert np.isclose(convergence, flux, atol=1e-9)
