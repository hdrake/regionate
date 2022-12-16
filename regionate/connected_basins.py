import numpy as np

from .utilities import *

class ConnectedBasins():
    def __init__(self, Basins):
        self.Basins = Basins
        self.Nb = len(self.Basins)
    
    def find_all_overlaps(self, tol=5.e3):
        self.overlaps = {}
        print("Finding overlaps for basin: ", end="")
        for i, (B1name, B1) in enumerate(self.Basins.items()):
            print(B1name, end=", ")
            for j, (B2name, B2) in enumerate(self.Basins.items()):
                if j<i:
                    overlaps = find_overlaps(B1, B2, tol=tol)
                    if len(overlaps[B1name]):
                        self.overlaps[(B1name, B2name)] = overlaps
                                
    def align_boundaries_with_overlap_sections(self):
        for Bname, B in self.Basins.items():
            overlap_list = [o for o in list(self.overlaps) if Bname in o]
            arbitrary_o = self.overlaps[overlap_list[0]][Bname]
            self.Basins[Bname] = roll_basin_boundary_to_align_with_overlap(B, arbitrary_o)
            
    def add_shared_corner_points(self):
        l = list(self.overlaps)
        corners = []

        for o in l:
            all_shared = [oo for oo in l if ((o[0] in oo) != (o[1] in oo))]
            for a in all_shared:
                if (((a[0], o[0]) in l) and ((a[0], o[1]) in l)):
                    corners.append([a[0], o[0], o[1]])
                if (((a[1], o[0]) in l) and ((a[1], o[1]) in l)):
                    corners.append([a[0], o[0], o[1]])

        corners = unique_list(corners)

        for c in corners:
            B = self.Basins[c[0]]

            idx01 = self.overlaps[(c[0], c[1])][c[0]]
            idx02 = self.overlaps[(c[0], c[2])][c[0]]

            dists = [
                (np.mod(idx01[0] -idx02[-1], B.N)),
                (np.mod(idx02[0] -idx01[-1], B.N)),
            ]
            inner = np.argmin(dists)

            if inner:
                coord = (B.lons[idx02[0]], B.lats[idx02[0]])
            else:
                coord = (B.lons[idx01[0]], B.lats[idx01[0]])

            # Find existing point closest to the corner and add it
            for cc in c:
                b = self.Basins[cc]
                if coord not in [(lon, lat) for (lon, lat) in zip(b.lons, b.lats)]:
                    dist = np.sqrt((b.lons - coord[0])**2 + (b.lats - coord[1])**2)
                    idx_insert = np.argmin(dist)
                    b.lons = np.concatenate((b.lons[:idx_insert], np.array([coord[0]]), b.lons[idx_insert:]))
                    b.lats = np.concatenate((b.lats[:idx_insert], np.array([coord[1]]), b.lats[idx_insert:]))
                    b.N = len(b.lons)
                    
def find_overlaps(B1, B2, fix_points_exact=False, tol=5.e3):
    overlaps = {B1.name:[], B2.name:[]}
    for i in range(B1.N):
        for j in range(B2.N):
            if haversine(B1.lons[i], B1.lats[i], B2.lons[j], B2.lats[j]) < tol:
                overlaps[B1.name].append(i)
                overlaps[B2.name].insert(0, j)
                if fix_points_exact:
                    B2.lons[j] = B1.lons[i]
                    B2.lats[j] = B1.lats[i]
    overlaps[B1.name] = consecutive_lists(overlaps[B1.name], mod=B1.N)
    overlaps[B2.name] = consecutive_lists(overlaps[B2.name], mod=B2.N)
    return group_overlaps(overlaps, B1, B2, tol=tol)


def group_overlaps(overlaps, B1, B2, tol=5.e3):
    grouped_overlaps = {B1.name: {}, B2.name: {}}
    group_num = 0
    for o2 in overlaps[B2.name]:
        for o1 in overlaps[B1.name]:
            if np.any(o1) and np.any(o2):
                for i in o1:
                    if np.any(haversine(B1.lons[i], B1.lats[i], B2.lons[o2], B2.lats[o2]) < tol):
                        grouped_overlaps[B1.name][group_num] = o1
                        grouped_overlaps[B2.name][group_num] = o2
                        group_num +=1
                        break
    return grouped_overlaps

def roll_basin_boundary_to_align_with_overlap(B, oidx):
    roll_idx = -oidx[0][0]
    new_oidx = {onum: np.mod(np.array(o) + roll_idx, B.N) for (onum, o) in oidx.items()}
    
    B.lons = np.roll(B.lons, roll_idx)
    B.lats = np.roll(B.lats, roll_idx)
    
    for onum, o in new_oidx.items():
        gaps = np.array([i for i in range(np.min(o), np.max(o)+1) if i not in o])
        if any(gaps):
            B.lons[gaps] = np.nan
            B.lats[gaps] = np.nan

    nan_idx = np.isnan(B.lons) | np.isnan(B.lats)
    B.lons = B.lons[~nan_idx] 
    B.lats = B.lats[~nan_idx]
    B.N = len(B.lons)
        
    return B