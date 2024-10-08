import numpy as np
from sectionate import distance_on_unit_sphere

def split_non_consecutive_list(data, mod=np.inf):
    """TO DO
    
    See also
    --------
    consecutive_lists, wrap_non_consecutive_listsoflists
    """
    if len(data) == 0:
        return [[]]
    else:
        data = iter(np.array(data).astype("int64"))
        val = next(data)
        chunk = []
        try:
            while True:
                chunk.append(val)
                val = next(data)
                if val != np.mod(chunk[-1] + 1, mod):
                    yield chunk
                    chunk = []
        except StopIteration:
            if chunk:
                yield chunk

def wrap_non_consecutive_listsoflists(lol, mod=np.inf):
    """TO DO"""
    if lol[0][0] == np.mod(lol[-1][-1]+1, mod):
        if len(lol)==1:
            lol[0] = lol[0][1:]
        else:
            lol = [lol[-1]+lol[0]] + lol[1:-1]
    lol = [l for l in lol if len(l)>1]
    return lol

def consecutive_lists(data, mod=np.inf):
    """TO DO"""
    if any(data):
        lol = [i for i in split_non_consecutive_list(data, mod=mod)]
        if any(lol):
            return wrap_non_consecutive_listsoflists(lol, mod=mod)
        else:
            return [[]]
    else:
        return [[]]
    
def coord_list(lons, lats):
    """Turns iterable longitudes and latitudes into a list of coordinate pairs"""
    return [(lon, lat) for (lon, lat) in zip(lons, lats)]
    
def unique_list(l):
    """Removes duplicate entries from list"""
    u = []
    for i, e in enumerate(l):
        # check if exists in unique_list or not
        if e not in u:
            u.append(e)
    return u

def unique_coords(coords, closeness_threshold=5.e3):
    """Removes duplicate coordinates from a list of (lon, lat) coordinate pairs

    PARAMETERS
    ----------
    coords : list of (lon, lat) coordinate pairs
    closeness_threshold : cutoff distance between "identical" points, as float (default: 5.e3)

    RETURNS
    -------
    coords

    See also
    --------
    unique_lonlat
    """
    lons = [coord[0] for coord in coords]
    lats = [coord[1] for coord in coords]
    uc = []
    for i, coord in enumerate(coords):
        # check if exists in unique_list or not
        if not(np.any( distance_on_unit_sphere(coord[0], coord[1], [c[0] for c in uc], [c[1] for c in uc]) < closeness_threshold )):
            uc.append(coord)
    return uc

def unique_lonlat(lons, lats, closeness_threshold=5.e3):
    """Removes duplicate coordinates from geographical coordinate arrays

    PARAMETERS
    ----------
    lons : array of longitudes in degrees
    lats : array of latitudes in degrees
    closeness_threshold : cutoff distance between "identical" points, as float (default: 5.e3)

    RETURNS
    -------
    lons, lats

    See also
    --------
    unique_coords
    """
    uc = unique_coords(coord_list(lons, lats), closeness_threshold=closeness_threshold)
    return np.array([c[0] for c in uc]), np.array([c[1] for c in uc])

def lon_mod(lon, lon_ref):
    """Longitudes modulo 360 degrees

    PARAMETERS
    ----------
    lon : input numpy array of longitudes [in degrees]
    lon_ref : reference longitude [in degrees]
    
    RETURNS
    -------
    np.ndarray with shape of `lon` but modulo 360 degrees
    
    """
    return lon + 360. *np.round((lon_ref - lon)/360.)

def loop(x):
    """Loop a list or numpy array by appending the first value to the end."""
    return np.append(x, x[0])

def flatten_dict_values(d):
    """
    Flattens iterable dictionary values into a single list.
    
    See also
    --------
    flatten_lol
    """
    return flatten_lol(d.values)

def flatten_lol(lol):
    """Flattens a list of lists into a single list"""
    return [i for l in lol for i in l]