# regionate

**Create xgcm-grid-consistent regional masks and boundaries for ocean model output.**

**regionate** turns a geographic polygon into a discrete, grid-aware region: it
snaps the polygon to a model grid and produces a boolean cell **mask** plus the
staggered `(u, v)` velocity **faces** that trace the region's boundary. Because
the mask and the boundary faces come from the same construction, a volume, mass,
or heat budget integrated over the masked region is exactly consistent with the
fluxes through its boundary — the property that makes closed regional budgets
possible.

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
contributing
api
```
