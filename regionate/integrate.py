import xarray as xr

def check_global_coverage(regions):
    total_mask = xr.zeros_like(list(regions.regions.values())[0].mask)
    for r in regions.regions.values():
        total_mask += r.mask
    if (total_mask == 1).sum() != total_mask.size:
        ValueError(f"Region {r.name} has incomplete or imperfect global coverage.")