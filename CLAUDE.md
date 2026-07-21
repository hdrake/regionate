# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`regionate` builds **xgcm-grid-consistent** regional masks and boundaries for ocean/climate model output. It supports **arbitrary `xgcm.Grid` topologies** — single-tile periodic and bipolar-fold (Arctic) MOM6 grids, and genuinely multi-tile grids defined by `face_connections` (e.g. ECCOv4r4 lat-lon-cap, cubed-sphere) — by driving all grid logic from the topology-aware `sectionate` API rather than hard-coded MOM6 specifics. Given a geographic polygon, it snaps the polygon to a discrete model grid, producing a boolean cell mask plus the staggered (u,v) velocity faces that trace the region's boundary — so that volume/mass/heat budgets integrated over the masked region are exactly consistent with fluxes through the boundary faces. It leans heavily on its sibling package [`sectionate`](https://github.com/MOM6-community/sectionate) for the section/face-tracing math (requires the topology-driven API on [`hdrake/sectionate@topology-driven-neighbors`](https://github.com/hdrake/sectionate/tree/topology-driven-neighbors), formerly PRs #47/#48), and on the bipolar north-fold boundary (`padding={..., "Y": {"fold": ...}}`, formerly xgcm#711) plus the multi-tile `face_connections` padding fix (xgcm#712) and the vector-pad fix (xgcm#749), all released in `xgcm >= 0.10.1` (on PyPI).

## AI Usage Policy

This repository has no formal AI Usage Policy of its own yet, so follow the one drafted for the sibling project xgcm ([`hdrake/xgcm@add-Claude.md`](https://github.com/hdrake/xgcm/tree/add-Claude.md)): **the person running the AI is responsible for every change it makes.** That means AI assistance must be disclosed, and the human must be able to explain the full diff. In practice:

- Do not produce changes the user could not stand behind and explain — no cargo-culted edits, no "it probably works."
- Surface uncertainty explicitly rather than hiding it behind confident-looking code.
- Keep diffs small and reviewable so they *can* be explained.

## Engineering norms

Adapted from the xgcm draft above; apply these when changing regionate's code.

1. **Deprecate by removing, not warning.** When you rename or remove public API, do not keep the old name limping along behind a `DeprecationWarning`. Remove it and make the old name fail immediately with a clear message (e.g. `raise ValueError("Argument 'old' renamed to 'new'.")`). Breaking changes are acceptable here; multi-release deprecation machinery is avoidable maintenance cost.
2. **Raise on bad input; never return a wrong answer.** regionate already does this — e.g. the corner-coordinate / grid validation the constructors run before building a mask. Keep it up: a silently wrong mask or boundary is worse than an exception in scientific software. On ambiguous/invalid/unsupported input, raise a specific error rather than guessing.
3. **Do not grow core dependencies.** The core deps are `sectionate` (the sibling that owns grid-section tracing and transports), plus `pyproj`, `geopandas`, `regionmask`, and `contourpy` (see `pyproject.toml`). Prefer pushing grid/transport logic *down into `sectionate`* over adding new dependencies here. Before adding any import, confirm it isn't pulling in a heavy new dependency.
4. **Every code change ships with tests.** Build test data from the synthetic grid fixtures in `regionate/tests/` (e.g. `initialize_spherical_grid()`) rather than hand-rolling datasets. For a bug fix, first write a test that fails, then make it pass. Exercise both construction directions (boundary→mask via `GriddedRegion`, mask→boundary via `MaskRegions`) and both orientations (clockwise/counterclockwise) where relevant.
5. **Version by effort; breaking changes are allowed.** Bump `regionate/version.py` (`__version__`, read by hatchling). Favor a clean break (norm 1) over a compatibility shim.

## Commands

```bash
pytest                                          # run the test suite
pytest regionate/tests/test_gridded_regions.py  # run the one test module
pytest -k <name>                                # run a single test by name

pip install -e .                                # editable install (after creating env)
```

Dev environment is conda-based (see README). CI (`.github/workflows/ci.yml`) installs `ci/environment.yml`, installs `sectionate` from `hdrake/sectionate@topology-driven-neighbors` (still unreleased), does `pip install -e .` (which pulls `xgcm >= 0.10.1` straight from PyPI — 0.10.1 ships the north-fold boundary and the multi-tile padding fixes), then runs `pytest` across Python 3.11–3.14. There is no linter configured. Optional end-to-end tests against real MOM6/ECCO grids are skipped unless `REGIONATE_REALDATA_TESTS=1`.

## Architecture

The class hierarchy in `regionate/region.py` is the spine of the package:

- **`Region`** — a named polygon, just `(lons, lats)` corner arrays. Can force counterclockwise winding (`is_section_counterclockwise` from sectionate) and prune duplicate points. No grid awareness.
- **`GriddedRegion(Region)`** — a `Region` bound to an `xgcm.Grid`. Construction calls `get_region_boundary_grid_indices` + `mask_from_grid_boundaries` (in `grid_conform.py`) to compute: corner grid indices (`i_c`, `j_c`, and the per-corner face index `f_c` — `None` for single-tile grids), the boundary `(u,v)` velocity faces, and the boolean cell `mask`. This is the central object — masks and boundaries are the whole point.
- **`BoundedRegion(GriddedRegion)`** — a GriddedRegion whose boundary is defined by named sections rather than a raw polygon.

Supporting modules:
- `regions.py` — **`Regions`**, a dict-like collection mapping names → `Region`/`GriddedRegion`; `overlaps.py` handles intersections between members.
- `grid_conform.py` — the core grid-snapping logic (polygon → grid indices → mask). `mask_from_grid_boundaries` rasterizes a boundary by splitting it at the antimeridian into a `[-180,180]` MultiPolygon and OR-ing per-piece `regionmask` masks (pole-encircling boundaries are extended to the South Pole first).
- `geometry.py` — antimeridian / MultiPolygon shapely helpers (`split_at_antimeridian`, `normalize_lon`, ...) used by `grid_conform.py`.
- `boundaries.py` — `grid_boundaries_from_mask` (inverse: recover boundary faces from a mask). A shared front-end `_trace_and_drop` traces the mask per face with `contourpy` and drops boundary segments interior to a seam (both separated cells in-mask in the topology-aware halo `_pad_center`); two back-ends then stitch the surviving arcs: `_single_tile_boundaries_from_mask` stitches by physical seam coincidence (handles walls, **periodic** axes, and the bipolar **fold** uniformly; `f_c=None`), and `_multitile_boundaries_from_mask` stitches on the grid's outer (shared-corner) corner-node graph from `sectionate.gridutils.outer_topology`: every traced corner resolves to a physical corner node (uniform across rotated/reversed seams, cube-vertex junctions, the pole, and grid cuts/folds), edges traced from both sides of an undeclared cut/fold annihilate as interior, surviving segments chain into closed loops by node identity, and each node is emitted as the native corner that stores it (`f_c_list`; `None` for single-tile). Boundaries through points stored on no face (LLC90's 4th Arctic-cap vertex, single-sided contact with its lon=-115 Antarctic cut) raise rather than fabricate indices. A region wrapping any seam yields one seam-consistent boundary (an annulus/strip may yield several loops); `_pad_center` pads only genuine seams (periodic/fold/face_connections) with real neighbours and everything else with NaN, so walls are never mistaken for seams. For the budget's *area side* on multi-tile staggered grids, use `outer_topology(grid).padded_transports(u, v)` — exact and xgcm-independent, and the only path that handles edges stored on no face (walls, grid cuts, a cap's un-stored vertex). (xgcm#749 fixes the bare-`DataArray` `grid.diff(..., other_component=)` path across rotated/reversed seams — the dict form was always exact — but even a correct pad cannot supply a halo for edges stored on no face, so `padded_transports` remains preferred here.)
- `integrate.py` — `check_global_coverage` validates that a `Regions` set is non-overlapping and tiles the globe.
- `utilities.py` — shared helpers (re-exported widely via `from .utilities import *`).

`regionate/__init__.py` flattens everything into the top-level namespace via `import *`, so public symbols are referenced as `regionate.Region`, `regionate.GriddedRegion`, etc.

### The `.gr` save format

`GriddedRegion.to_gr(path)` writes a **directory** (`<name>.gr/`), not a single file: `region.nc` (boundary coords `lons_c`/`lats_c`, corner indices `i_c`/`j_c`/`f_c`, optional `lons_uv`/`lats_uv`, and the `mask`) plus the grid dataset and child-section sub-directories as separate NetCDFs. `open_gr(path, ds_to_grid)` reloads it — note it takes a `ds_to_grid` callback because the `xgcm.Grid` cannot be serialized directly and must be reconstructed from the saved dataset. The face index `f_c` is persisted and reloaded so multi-tile regions do not silently reload as single-tile.

## Conventions

- Suffixes encode grid staggering: `_c` = tracer/corner points, `_uv` = velocity faces. Preserve this when adding coordinates.
- Many functions take an `xgcm.Grid` instance as `grid` and assume `grid._ds` holds the underlying xarray dataset with MOM6-style symmetric staggering.
- `examples/` contains runnable Jupyter notebooks (thickness/heat/mass budgets) that double as the worked-example documentation; `examples/load_example_model_grid*.py` build the demo grids from the NetCDFs in `data/`.
