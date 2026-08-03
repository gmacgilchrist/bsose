import xarray as xr
import os

basedir = '/gws/ssde/j25a/co2clim/datasets/bSOSE/ITER139'
outdir = '/gws/ssde/j25a/co2clim/public/graemem/datasets/bSOSE/subset_WeddellGyre'
filename_suffix = '_bsoseI139_2013to2021_5dy.nc'
variables = ['MLD','Theta','Salt','UVel','VVel','WVel']

for variable in variables:
    print(variable)
    filename = variable+filename_suffix
    inpath = '/'.join([basedir,filename])
    outpath = '/'.join([outdir,filename])
    # open dataset
    ds = xr.open_dataset(inpath)
    # wrap x-coordinate
    ds_rolled = ds.roll(XC=round(len(ds["XC"])/2), roll_coords=True)
    shift = ds_rolled["XC"][0]
    tmp = (ds_rolled["XC"]+180)%360-shift
    ds_rolled["XC"] = tmp
    # clip selection
    select = {"XC":slice(-70,110),"YC":slice(-90,-50)}
    iselect = {'time':slice(-2*73,None)}
    ds_clipped = ds_rolled.sel(select).isel(iselect)
    # save to netcdf
    if os.path.exists(outpath):
        continue
    else:
            ds_clipped.to_netcdf(outpath)