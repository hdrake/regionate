import numpy as np
from .utilities import *

def find_indices_of_overlaps(b1, b2, closeness_threshold=5.e3, face_indices=False):
    overlaps = haversine(
            b1.lons[:, np.newaxis], b1.lats[:, np.newaxis],
            b2.lons[np.newaxis, :], b2.lats[np.newaxis, :]
        ) < closeness_threshold

    b1_oidx = consecutive_lists( np.where(np.any(overlaps, axis=1))[0], mod=b1.lons.size)
    b2_oidx = consecutive_lists( np.where(np.any(overlaps, axis=0))[0], mod=b2.lons.size)
    
    if face_indices:
        b1_oidx = [l[:-1] for l in b1_oidx]
        b2_oidx = [l[:-1] for l in b2_oidx]
    
    return group_overlaps({b1.name: b1_oidx, b2.name: b2_oidx}, b1, b2, closeness_threshold=closeness_threshold)

def group_overlaps(overlaps, b1, b2, closeness_threshold=5.e3):
    grouped_overlaps = {b1.name: {}, b2.name: {}}
    group_num = 0
    for o2 in overlaps[b2.name]:
        for o1 in overlaps[b1.name]:
            if np.any(o1) and np.any(o2):
                if np.any(haversine(
                        b1.lons[o1, np.newaxis], b1.lats[o1, np.newaxis],
                        b2.lons[np.newaxis, o2], b2.lats[np.newaxis, o2]
                    ) < closeness_threshold):
                    grouped_overlaps[b1.name][group_num] = o1
                    grouped_overlaps[b2.name][group_num] = o2
                    group_num +=1
    return grouped_overlaps


def align_boundaries_with_overlap_sections(region, remove_gaps=True):
    for bname, b in region.Basins.items():
        overlap_list = [o for o in list(region.overlaps) if bname in o]
        if len(overlap_list) != 0:
            arbitrary_o = region.overlaps[overlap_list[0]][bname]
            region.Basins[bname] = roll_boundary_to_align_with_overlap(
                b, arbitrary_o, remove_gaps=remove_gaps
            )
            
def roll_boundary_to_align_with_overlap(b, oidx, remove_gaps=True):
    roll_idx = -consecutive_lists(oidx[0], mod=b.lons.size)[0][0]
    new_oidx = {onum: np.mod(np.array(o) + roll_idx, b.lons.size) for (onum, o) in oidx.items()}
    
    b.lons = np.roll(b.lons, roll_idx)
    b.lats = np.roll(b.lats, roll_idx)
    
    if remove_gaps:
        for onum, o in new_oidx.items():
            gaps = np.array([i for i in range(np.min(o), np.max(o)+1) if i not in o])
            if any(gaps):
                b.lons[gaps] = np.nan
                b.lats[gaps] = np.nan

    nan_idx = np.isnan(b.lons) | np.isnan(b.lats)
    b.lons = b.lons[~nan_idx] 
    b.lats = b.lats[~nan_idx]
        
    return b