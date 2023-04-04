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
            self.counterclockwise = self.is_boundary_counterclockwise()
        if force_ccw:
            self.make_counterclockwise()

    def copy(self, remove_duplicate_points=False):
        return Basin(self.lons.copy(), self.lats.copy(), self.name, counterclockwise=self.counterclockwise, remove_duplicate_points=remove_duplicate_points)
    
    def is_boundary_counterclockwise(self):
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