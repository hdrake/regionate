# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`regionate` builds **xgcm-grid-consistent** regional masks and boundaries for ocean/climate model output. It supports **arbitrary `xgcm.Grid` topologies** — single-tile periodic and bipolar-fold (Arctic) MOM6 grids, and genuinely multi-tile grids defined by `face_connections` (e.g. ECCOv4r4 lat-lon-cap, cubed-sphere) — by driving all grid logic from the topology-aware `sectionate` API rather than hard-coded MOM6 specifics. Given a geographic polygon, it snaps the polygon to a discrete model grid, producing a boolean cell mask plus the staggered (u,v) velocity faces that trace the region's boundary — so that volume/mass/heat budgets integrated over the masked region are exactly consistent with fluxes through the boundary faces. It leans heavily on its sibling package [`sectionate`](https://github.com/MOM6-community/sectionate) for the section/face-tracing math (requires the topology-driven [sectionate#47](https://github.com/MOM6-community/sectionate/pull/47)).

## Commands

```bash
pytest                                          # run the test suite
pytest regionate/tests/test_gridded_regions.py  # run the one test module
pytest -k <name>                                # run a single test by name

pip install -e .                                # editable install (after creating env)
```

Dev environment is conda-based (see README). CI (`.github/workflows/ci.yml`) installs `ci/environment.yml`, installs `sectionate` from the topology-driven PR ref (`refs/pull/47/head`) until it is released, does `pip install -e .`, then runs `pytest` across Python 3.11–3.14. There is no linter configured. Optional end-to-end tests against real MOM6/ECCO grids are skipped unless `REGIONATE_REALDATA_TESTS=1`.

## Architecture

The class hierarchy in `regionate/region.py` is the spine of the package:

- **`Region`** — a named polygon, just `(lons, lats)` corner arrays. Can force counterclockwise winding (`is_section_counterclockwise` from sectionate) and prune duplicate points. No grid awareness.
- **`GriddedRegion(Region)`** — a `Region` bound to an `xgcm.Grid`. Construction calls `get_region_boundary_grid_indices` + `mask_from_grid_boundaries` (in `grid_conform.py`) to compute: corner grid indices (`i_c`, `j_c`, and the per-corner face index `f_c` — `None` for single-tile grids), the boundary `(u,v)` velocity faces, and the boolean cell `mask`. This is the central object — masks and boundaries are the whole point.
- **`BoundedRegion(GriddedRegion)`** — a GriddedRegion whose boundary is defined by named sections rather than a raw polygon.

Supporting modules:
- `regions.py` — **`Regions`**, a dict-like collection mapping names → `Region`/`GriddedRegion`; `overlaps.py` handles intersections between members.
- `grid_conform.py` — the core grid-snapping logic (polygon → grid indices → mask). `mask_from_grid_boundaries` rasterizes a boundary by splitting it at the antimeridian into a `[-180,180]` MultiPolygon and OR-ing per-piece `regionmask` masks (pole-encircling boundaries are extended to the South Pole first).
- `geometry.py` — antimeridian / MultiPolygon shapely helpers (`split_at_antimeridian`, `normalize_lon`, ...) used by `grid_conform.py`.
- `boundaries.py` — `grid_boundaries_from_mask` (inverse: recover boundary faces from a mask). Single-tile grids trace with `contourpy`; multi-tile grids trace per-face and stitch across tile seams using the grid's `face_connections` topology, returning a per-corner face index (`f_c_list`).
- `integrate.py` — `check_global_coverage` validates that a `Regions` set is non-overlapping and tiles the globe.
- `utilities.py` — shared helpers (re-exported widely via `from .utilities import *`).

`regionate/__init__.py` flattens everything into the top-level namespace via `import *`, so public symbols are referenced as `regionate.Region`, `regionate.GriddedRegion`, etc.

### The `.gr` save format

`GriddedRegion.to_gr(path)` writes a **directory** (`<name>.gr/`), not a single file: `region.nc` (boundary coords `lons_c`/`lats_c`, corner indices `i_c`/`j_c`/`f_c`, optional `lons_uv`/`lats_uv`, and the `mask`) plus the grid dataset and child-section sub-directories as separate NetCDFs. `open_gr(path, ds_to_grid)` reloads it — note it takes a `ds_to_grid` callback because the `xgcm.Grid` cannot be serialized directly and must be reconstructed from the saved dataset. The face index `f_c` is persisted and reloaded so multi-tile regions do not silently reload as single-tile.

## Conventions

- Suffixes encode grid staggering: `_c` = tracer/corner points, `_uv` = velocity faces. Preserve this when adding coordinates.
- Many functions take an `xgcm.Grid` instance as `grid` and assume `grid._ds` holds the underlying xarray dataset with MOM6-style symmetric staggering.
- `examples/` contains runnable Jupyter notebooks (thickness/heat/mass budgets) that double as the worked-example documentation; `examples/load_example_model_grid*.py` build the demo grids from the NetCDFs in `data/`.
