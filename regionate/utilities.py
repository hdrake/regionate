import numpy as np

def haversine(lon1, lat1, lon2, lat2, R=6.371e6):
    lon1, lat1, lon2, lat2 = np.deg2rad(lon1), np.deg2rad(lat1), np.deg2rad(lon2), np.deg2rad(lat2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))

def unique_list(l):
    u = []
    for i, e in enumerate(l):
        # check if exists in unique_list or not
        if e not in u:
            u.append(e)
    return u

def split_non_consecutive_list(data, mod=np.inf):
    data = iter(np.array(data).astype("int64"))
    val = next(data)
    chunk = []
    try:
        while True:
            chunk.append(val)
            val = next(data)
            if ((val != np.mod(chunk[-1] + 1, mod)) and (val != np.mod(chunk[-1] + 2, mod)) and (val != np.mod(chunk[-1] + 3, mod))):
                yield chunk
                chunk = []
    except StopIteration:
        if chunk:
            yield chunk

def wrap_non_consecutive_listsoflists(lol, mod=np.inf):
    if lol[0][0] == np.mod(lol[-1][-1]+1, mod):
        if len(lol)==1:
            lol[0] = lol[0][1:]
        else:
            lol = [lol[-1]+lol[0]] + lol[1:-1]
    lol = [l for l in lol if len(l)>1]
    return lol

def consecutive_lists(data, mod=np.inf):
    if any(data):
        lol = [i for i in split_non_consecutive_list(data, mod=mod)]
        if any(lol):
            return wrap_non_consecutive_listsoflists(lol, mod=mod)
        else:
            return [[]]
    else:
        return [[]]