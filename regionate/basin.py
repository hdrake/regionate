import numpy as np

from .utilities import *

class Basin:
    def __init__(self, lons, lats, name, counterclockwise=None, force_ccw=False, remove_duplicate_points=False):
        self.name = name
        self.lons = lons
        self.lats = lats
        
        if remove_duplicate_points:
            self.remove_duplicate_points()

        if counterclockwise is not None:
            self.counterclockwise = counterclockwise
        else:
            self.counterclockwise = self.check_if_boundary_clockwise()
        if force_ccw:
            self.make_counterclockwise()

    def copy(self, remove_duplicate_points=False):
        return Basin(self.lons.copy(), self.lats.copy(), self.name, counterclockwise=self.counterclockwise, remove_duplicate_points=remove_duplicate_points)
    
    def check_if_boundary_clockwise(self):
        
        circumgeo = (np.abs(np.sum(np.diff(self.lons))) >= 180.) | (np.abs(np.sum(np.diff(self.lats))) >= 180.)
        if circumgeo:
            raise ValueError("""Keyword argument counterclockwise must be explicitly set for sections
            that circumnavigate the globe, because orientations are otherwise ambiguous.""")
        else:
            lons = np.append(self.lons, self.lons[0])
            lats = np.append(self.lats, self.lats[0])
            signed_area = 0.
            for i in range(self.lons.size):
                signed_area += (lons[i+1]-lons[i])*(lats[i+1]+lats[i])
            return signed_area < 0.
    
    def make_counterclockwise(self):
        if not(self.counterclockwise):
            self.lons = self.lons[::-1]
            self.lats = self.lats[::-1]
            self.counterclockwise = True

    def remove_duplicate_points(self, closeness_threshold=5.e3):
        self.lons, self.lats = unique_lonlat(self.lons, self.lats, closeness_threshold=closeness_threshold)
        
def prune_redundant_boundary_points(b, closeness_threshold=5.e3):
    """ Find repeated points, delete all indices in shortest path between them,
    and remove the second visit to the point. Update valid sectionate points.""" 
    repeats = haversine(
        b.lons[:, np.newaxis], b.lats[:, np.newaxis],
        b.lons[np.newaxis, :], b.lats[np.newaxis, :]
    ) < closeness_threshold

    duals = []
    repeat_idx, repeat_idx_dual = np.where( repeats - np.identity(b.lons.size) )

    keep_flags = np.ones(b.lons.shape, dtype=np.bool_)
    for n, (i1, i2) in enumerate(zip(repeat_idx, repeat_idx_dual)):
        if i1 in duals: continue
        duals.append(i2)
        if np.abs(i2 - i1) < b.lons.size//2:
            keep_idx = np.arange(i1, i2)
        else:
            keep_idx = np.append(np.arange(0, i1), np.arange(i2, b.lons.size))
        keep_flags[keep_idx] = 0.
        
    b.lons = b.lons[keep_flags]
    b.lats = b.lats[keep_flags]
    b.lons_uv = b.lons_uv[keep_flags]
    b.lats_uv = b.lats_uv[keep_flags]
    b.i = np.append(b.i[:-1][keep_flags], b.i[:-1][keep_flags][0])
    b.j = np.append(b.j[:-1][keep_flags], b.j[:-1][keep_flags][0])