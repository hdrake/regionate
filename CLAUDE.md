# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`regionate` builds **xgcm-grid-consistent** regional masks and boundaries for ocean/climate model output (currently mostly MOM6, but with plans to support any ocean model). Given a geographic polygon, it snaps the polygon to a discrete model grid, producing a boolean cell mask plus the staggered (u,v) velocity faces that trace the region's boundary — so that volume/mass/heat budgets integrated over the masked region are exactly consistent with fluxes through the boundary faces. It leans heavily on its sibling package [`sectionate`](https://github.com/MOM6-community/sectionate) for the section/face-tracing math.

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
5. **Version by effort; breaking changes are allowed.** The version comes from the git tag, not from a file in the tree — see "Versioning" below, and do not add a literal back to `regionate/version.py`. Favor a clean break (norm 1) over a compatibility shim.

## Commands

```bash
pytest                                          # run the test suite
pytest regionate/tests/test_gridded_regions.py  # run the one test module
pytest -k <name>                                # run a single test by name

pip install -e .                                # editable install (after creating env)
```

Dev environment is conda-based (see README). CI (`.github/workflows/ci.yml`) installs `ci/environment.yml`, does `pip install -e .`, then runs `pytest` across Python 3.11–3.14. There is no linter configured.

## Versioning

**The git tag is the single source of truth.** `hatch-vcs` (`[tool.hatch.version] source = "vcs"`) derives the version from the tag at build time and writes it to `regionate/_version.py`, which is **gitignored** — there is no version string in the source tree. `regionate/version.py` is a thin shim that imports from it, with a `0.0.0+unknown` fallback for an un-built checkout.

Consequences worth remembering when editing:

- **Never add a version literal back to the tree**, and never "fix" a `0.0.0+unknown` by hardcoding one — it means the package was imported without being built or installed.
- **Any CI job that installs the package needs `fetch-depth: 0`.** A shallow clone cannot see the tag, so hatch-vcs silently resolves a `0.1.devN` version instead of failing. The checkouts in `ci.yml` and `publish-to-pypi.yml` set it, and `.readthedocs.yaml` unshallows in `post_checkout` for the same reason. The same applies to a `pip install git+https://…` of a fork with no tags: it reports `0.1.devN`, which can fall below a downstream floor.
- `_version.py` **is** shipped inside the sdist, so building from the sdist (as conda-forge does) works with no git present. Do not add it to `[tool.hatch.build] exclude`.
- conda builds run with `--no-build-isolation`, so `hatch-vcs` must sit next to `hatchling` in the recipe's `host` requirements — in `conda/meta.yaml` here and on the conda-forge feedstock.
- Releasing is just publishing a GitHub Release tagged `vX.Y.Z`; there is no bump commit. See "Releasing" in `README.md`.

## Architecture

The class hierarchy in `regionate/region.py` is the spine of the package:

- **`Region`** — a named polygon, just `(lons, lats)` corner arrays. Can force counterclockwise winding (`is_section_counterclockwise` from sectionate) and prune duplicate points. No grid awareness.
- **`GriddedRegion(Region)`** — a `Region` bound to an `xgcm.Grid`. Construction calls `get_region_boundary_grid_indices` + `mask_from_grid_boundaries` (in `grid_conform.py`) to compute: corner grid indices (`i_c`, `j_c`), the boundary `(u,v)` velocity faces, and the boolean cell `mask`. This is the central object — masks and boundaries are the whole point.
- **`BoundedRegion(GriddedRegion)`** — a GriddedRegion whose boundary is defined by named sections rather than a raw polygon.

Supporting modules:
- `regions.py` — **`Regions`**, a dict-like collection mapping names → `Region`/`GriddedRegion`; `overlaps.py` handles intersections between members.
- `grid_conform.py` — the core grid-snapping logic (polygon → grid indices → mask), using `regionmask`, `shapely`, and sectionate.
- `boundaries.py` — `grid_boundaries_from_mask` (inverse: recover boundary faces from a mask).
- `integrate.py` — `check_global_coverage` validates that a `Regions` set is non-overlapping and tiles the globe.
- `utilities.py` — shared helpers (re-exported widely via `from .utilities import *`).

`regionate/__init__.py` flattens everything into the top-level namespace via `import *`, so public symbols are referenced as `regionate.Region`, `regionate.GriddedRegion`, etc.

### The `.gr` save format

`GriddedRegion.save()` writes a **directory** (`<name>.gr/`), not a single file: `region.nc` plus the grid dataset and section sub-directories as separate NetCDFs. `open_gr(path, ds_to_grid)` reloads it — note it takes a `ds_to_grid` callback because the `xgcm.Grid` cannot be serialized directly and must be reconstructed from the saved dataset.

## Conventions

- Suffixes encode grid staggering: `_c` = tracer/corner points, `_uv` = velocity faces. Preserve this when adding coordinates.
- Many functions take an `xgcm.Grid` instance as `grid` and assume `grid._ds` holds the underlying xarray dataset with MOM6-style symmetric staggering.
- `examples/` contains runnable Jupyter notebooks (thickness/heat/mass budgets) that double as the worked-example documentation; `examples/load_example_model_grid*.py` build the demo grids from the NetCDFs in `data/`.
