import numpy as np

from .utilities import *

class Basin:
    def __init__(self, lons, lats, name, force_ccw=True):
        self.name = name
        self.lons = lons
        self.lats = lats
        self.remove_duplicate_points()

        self.section_is_clockwise = self.check_if_boundary_clockwise()
        self.make_counterclockwise()
    
    def check_if_boundary_clockwise(self):
        lons = np.append(self.lons, self.lons[0])
        lats = np.append(self.lats, self.lats[0])
        signed_area = 0.
        for i in range(self.N):
            signed_area += (lons[i+1]-lons[i])*(lats[i+1]+lats[i])
        return signed_area >= 0.
    
    def make_counterclockwise(self):
        if self.section_is_clockwise:
            self.lons = self.lons[::-1]
            self.lats = self.lats[::-1]
            
    def remove_duplicate_points(self, tol=5.e3):
        coords = [(lon, lat) for (lon, lat) in zip(self.lons, self.lats)]

        unique_coords = []
        for i, c in enumerate(coords):
            # check if exists in unique_list or not
            unique_lons = np.array([lon for (lon, lat) in unique_coords])
            unique_lats = np.array([lat for (lon, lat) in unique_coords])
            if np.any( haversine(c[0], c[1], unique_lons, unique_lats) < tol ):
                self.lons[i] = np.nan
                self.lats[i] = np.nan
            else:
                unique_coords.append(c)
        nan_idx = np.isnan(self.lons) | np.isnan(self.lats)
        self.lons = self.lons[~nan_idx]
        self.lats = self.lats[~nan_idx]
        self.N = len(self.lons)