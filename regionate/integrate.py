import xarray as xr

def check_global_coverage(r):
    total_mask = xr.zeros_like(list(r.Basins.values())[0].mask)
    for b in r.Basins.values():
        total_mask += b.mask
    if (total_mask == 1).sum() != total_mask.size:
        ValueError(f"Region {r.name} has incomplete or imperfect global coverage.")