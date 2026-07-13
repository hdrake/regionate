import urllib.request
import shutil
import os
import xarray as xr
import xgcm

def download_MOM6_example_data(file_name):
    # download the data
    url = 'https://zenodo.org/record/15384910/files/'
    destination_path = f"../data/{file_name}"
    if not os.path.exists(destination_path):
        print(f"File '{file_name}' being downloaded to {destination_path}.")
        with urllib.request.urlopen(url + file_name) as response, open(destination_path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        print(f"File '{file_name}' has completed download to {destination_path}.")
    else:
        print(f"File '{file_name}' already exists at {destination_path}. Skipping download.")

    return destination_path

def load_MOM6_example_grid(file_name, fold=False):
    destination_path = download_MOM6_example_data(file_name)
    ds = xr.open_dataset(destination_path).fillna(0.)
    if "z_l" not in ds.dims:
        ds = ds.expand_dims(["z_l"]).assign_coords({
            "z_l":xr.DataArray([3000], dims=("z_l",)),
            "z_i":xr.DataArray([0,6000], dims=("z_i",))
        })
    return construct_grid(ds, fold=fold)

def load_MOM6_zint_mass_budget(fold=False):
    file_name = 'MOM6_global_example_vertically_integrated_mass_budget_v0_0_6.nc'
    return load_MOM6_example_grid(file_name, fold=fold)

def load_MOM6_zint_heat_budget(fold=False):
    file_name = 'MOM6_global_example_vertically_integrated_heat_budget_v0_0_6.nc'
    return load_MOM6_example_grid(file_name, fold=fold)

def construct_grid(ds, fold=False):
    coords={
        'X': {'center': 'xh', 'outer': 'xq'},
        'Y': {'center': 'yh', 'outer': 'yq'},
    }
    # This is a tripolar grid: its northern edge is a bipolar fold, not a wall. Pass
    # `fold=True` (requires an xgcm with north-fold support, hdrake/xgcm@dev-v1.0.0) to
    # declare it, so `regionate` traces regions straddling the Arctic fold into a single
    # boundary loop (see notebook 3). The default `Y='extend'` treats the fold as a wall,
    # which is adequate for regions away from the Arctic and splits fold-straddling ones in two.
    padding = {'X':'periodic', 'Y':({'fold':'corner'} if fold else 'extend')}
    metrics = {('X','Y'):'areacello'}
    grid = xgcm.Grid(ds, coords=coords, metrics=metrics, padding=padding, autoparse_metadata=False)
    return grid