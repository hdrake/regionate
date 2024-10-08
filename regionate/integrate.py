import xarray as xr

def check_global_coverage(regions):
    """Check whether the masks in the `Regions` instance are non-overlapping and provide complete global coverage"""
    total_mask = xr.zeros_like(list(regions.region_dict.values())[0].mask)
    for r in regions.region_dict.values():
        total_mask += r.mask
    if (total_mask == 1).sum() != total_mask.size:
        ValueError(f"Region {r.name} has incomplete or imperfect global coverage.")