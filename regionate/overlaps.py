import numpy as np
from .utilities import *

def find_indices_of_overlaps(r1, r2, closeness_threshold=5.e3, face_indices=False):
    overlaps = haversine(
            r1.lons[:, np.newaxis], r1.lats[:, np.newaxis],
            r2.lons[np.newaxis, :], r2.lats[np.newaxis, :]
        ) < closeness_threshold

    r1_oidx = consecutive_lists( np.where(np.any(overlaps, axis=1))[0], mod=r1.lons.size)
    r2_oidx = consecutive_lists( np.where(np.any(overlaps, axis=0))[0], mod=r2.lons.size)
    
    if face_indices:
        r1_oidx = [l[:-1] for l in r1_oidx]
        r2_oidx = [l[:-1] for l in r2_oidx]
    
    return group_overlaps({r1.name: r1_oidx, r2.name: r2_oidx}, r1, r2, closeness_threshold=closeness_threshold)

def group_overlaps(overlaps, r1, r2, closeness_threshold=5.e3):
    grouped_overlaps = {r1.name: {}, r2.name: {}}
    group_num = 0
    for o2 in overlaps[r2.name]:
        for o1 in overlaps[r1.name]:
            if np.any(o1) and np.any(o2):
                if np.any(haversine(
                        r1.lons[o1, np.newaxis], r1.lats[o1, np.newaxis],
                        r2.lons[np.newaxis, o2], r2.lats[np.newaxis, o2]
                    ) < closeness_threshold):
                    grouped_overlaps[r1.name][group_num] = o1
                    grouped_overlaps[r2.name][group_num] = o2
                    group_num +=1
    return grouped_overlaps


def align_boundaries_with_overlap_sections(regions, remove_gaps=True):
    for rname, r in regions.region_dict.items():
        overlap_list = [o for o in list(regions.overlaps) if rname in o]
        if len(overlap_list) != 0:
            arbitrary_o = regions.overlaps[overlap_list[0]][rname]
            regions.region_dict[rname] = roll_boundary_to_align_with_overlap(
                r, arbitrary_o, remove_gaps=remove_gaps
            )
            
def roll_boundary_to_align_with_overlap(r, oidx, remove_gaps=True):
    roll_idx = -consecutive_lists(oidx[0], mod=r.lons.size)[0][0]
    new_oidx = {onum: np.mod(np.array(o) + roll_idx, r.lons.size) for (onum, o) in oidx.items()}
    
    r.lons = np.roll(r.lons, roll_idx)
    r.lats = np.roll(r.lats, roll_idx)
    
    if remove_gaps:
        for onum, o in new_oidx.items():
            gaps = np.array([i for i in range(np.min(o), np.max(o)+1) if i not in o])
            if any(gaps):
                r.lons[gaps] = np.nan
                r.lats[gaps] = np.nan

    nan_idx = np.isnan(r.lons) | np.isnan(r.lats)
    r.lons = r.lons[~nan_idx] 
    r.lats = r.lats[~nan_idx]
        
    return r