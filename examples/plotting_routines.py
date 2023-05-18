import matplotlib.pyplot as plt

def geoplot(ds, da, *args, **kwargs):
    pc = plt.pcolormesh(
        ds['geolon_c'],
        ds['geolat_c'],
        da,
        *args,
        **kwargs
    )
    return pc