import numpy as np

from .basin import *
from .overlaps import *
from .utilities import *

class Region():
    def __init__(self, Basins, name=None):
        self.Basins = Basins
        if name is not None:
            self.name = name
    
    def find_all_overlaps(self, closeness_threshold=5.e3, face_indices=False):
        self.overlaps = {}
        for i, (b1name, b1) in enumerate(self.Basins.items()):
            for j, (b2name, b2) in enumerate(self.Basins.items()):
                if b1name<b2name:
                    overlaps = find_indices_of_overlaps(b1, b2, closeness_threshold=closeness_threshold, face_indices=face_indices)
                    if len(overlaps[b1name]):
                        self.overlaps[sorted_tuple((b1name, b2name))] = overlaps
                
    def copy(self, remove_duplicate_points=False):
        return Region({b.name: b.copy(remove_duplicate_points=remove_duplicate_points) for b in self.Basins.values()})

def sorted_tuple(s):
    return tuple(sorted(s, key=int))