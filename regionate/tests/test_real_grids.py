"""Optional end-to-end checks against real ocean-model grids.

These are skipped by default. Enable them by setting REGIONATE_REALDATA_TESTS=1
(and, for the ECCO check, pointing REGIONATE_ECCO_GRID at a NetCDF file holding
an ECCOv4r4 LLC90 grid with `face_connections`). They validate the multi-model
goal of the overhaul on genuine MOM6 (tripolar) and ECCO (lat-lon-cap) output.
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
    not os.path.isfile(ECCO_FILE),
    reason="ECCO LLC90 geometry not downloaded (see examples/load_example_ECCO_grid.py)",
)
def _load_ecco():
    import sys
    sys.path.insert(0, os.path.abspath(EXAMPLES_DIR))
    from load_example_ECCO_grid import load_ECCO_LLC90_grid, atlantic_basin_mask
    grid = load_ECCO_LLC90_grid(data_dir=os.path.abspath(DATA_DIR))
    return grid, atlantic_basin_mask


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


def test_ecco_atlantic_basin_obeys_discrete_divergence_theorem():
    """The whole point of the package: a region's budget must close against the
    fluxes through its traced boundary. On the real LLC90 grid, for the full
    Atlantic basin (spanning rotated seams), the net flux through every boundary
    velocity face (summed over all loops) equals the flux convergence summed over
    the masked cells -- to machine precision, for an arbitrary transport field.

    The cell-centred convergence is taken with xgcm's vector-aware
    ``grid.diff_2d_vector``: across a 90-degree LLC seam the U-component rotates into
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

    div = grid.diff_2d_vector({"X": umo, "Y": vmo}, to="center",
                              boundary="fill", fill_value=np.nan)
    convergence = float((-(div["X"] + div["Y"])).where(mask, 0.).sum())

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
