# regionate

**Create xgcm-grid-consistent regional masks and boundaries for ocean model output.**

**regionate** turns a geographic polygon into a discrete, grid-aware region: it
snaps the polygon to a model grid and produces a boolean cell **mask** plus the
staggered `(u, v)` velocity **faces** that trace the region's boundary. Because
the mask and the boundary faces come from the same construction, a volume, mass,
or heat budget integrated over the masked region is exactly consistent with the
fluxes through its boundary — the property that makes closed regional budgets
possible.

All grid logic is driven by the grid's own topology, so this works on **arbitrary
`xgcm.Grid` topologies**: single-tile periodic grids, bipolar-fold (Arctic) MOM6
grids, and genuinely multi-tile grids defined by `face_connections` (e.g. the
ECCOv4r4 lat-lon-cap grid, cubed spheres).

It leans on its sibling package
[`sectionate`](https://github.com/MOM6-community/sectionate) for the
section- and face-tracing math, and builds on
[`xgcm`](https://xgcm.readthedocs.io/en/stable/) grids.

The central object is `GriddedRegion` (a polygon bound to an `xgcm.Grid`);
`BoundedRegion` defines a region from named sections instead of a raw polygon,
and `Regions` collects many of them. See the worked examples below for thickness,
heat, and mass budgets over named regions.

```{toctree}
:maxdepth: 2
:caption: Contents

installation
examples/1_thickness_budget
examples/2_advective_heat_convergence
examples/3_Arctic_heat_CM4p25
examples/4_bounded_by_named_sections
examples/5_ECCO_LLC90_multiface_regions
examples/6_idealized_corner_cases
contributing
api
```
