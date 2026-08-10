import xarray as xr

def check_global_coverage(regions):
    """Raise if the region masks do not partition the domain.

    A valid partition covers every cell exactly once: each cell must be claimed by
    exactly one region's mask, with no gaps (cells claimed by none) and no overlaps
    (cells claimed by more than one).
    """
    masks = [r.mask for r in regions.region_dict.values()]
    coverage = xr.zeros_like(masks[0], dtype=int)
    for m in masks:
        coverage += m.astype(int)
    n_gap = int((coverage == 0).sum())
    n_overlap = int((coverage > 1).sum())
    if n_gap or n_overlap:
        raise ValueError(
            f"Region masks do not partition the domain: {n_gap} cell(s) covered by "
            f"no region and {n_overlap} cell(s) covered by more than one region."
        )