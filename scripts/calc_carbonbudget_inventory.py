import xarray as xr
import pandas as pd
from xgcm import Grid
import numpy as np
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar
import os

import bsose.preprocess as pp

ds,xgrid = pp.load_bsose()

# Define time metric
# HACK: trouble with time difference metric, so here just setting up own array with 5-days in seconds 
dt = xr.DataArray(432000*np.ones(shape=(438)),dims='time')
# Reference density
rho0 = 1035.0

# Define some masks
# Mask to set surface velocity point to zero
tmp = np.ones(len(ds['Zl']))
tmp[0]=0
maskZl = xr.DataArray(tmp,dims=['Zl'],coords={'Zl':ds['Zl']})

# Mask to set surface tracer point to one (everything else to zero)
tmp = np.ones(len(ds['Z']))
tmp[1:]=0
maskZ = xr.DataArray(tmp,dims=['Z'],coords={'Z':ds['Z']})

budget = xr.Dataset()
# Thicknesses
h = ds['drC']+ds['ETAN']*maskZ
eta = ds['ETAN']*maskZ

# TERMS
# Tendency
h_snaps = ds['ETAN_snaps']*maskZ + ds['drC']
hPhi = ds['TRAC01_snaps']*h_snaps
budget['TEND'] = xgrid.diff(hPhi,'T')/dt
# Advection
ADVc = -(xgrid.diff(ds['ADVxTr01'],'X')+
         xgrid.diff(ds['ADVyTr01'],'Y',boundary='extend')+
         (-1*xgrid.diff(ds['ADVrTr01'],'Z',boundary='extend')))/ds['vC']
budget['ADV'] = h*ADVc
# Diffusion
DIFFc = -(xgrid.diff(ds['DFxETr01'],'X')+
          xgrid.diff(ds['DFyETr01'],'Y',boundary='extend')+
          (-1*xgrid.diff(ds['DFrITr01'],'Z',boundary='extend')))/ds['vC']
budget['DIFF'] = h*DIFFc
# Air-sea flux
SURFc = maskZ*(ds['BLGCFLX']/ds['drC'])
budget['SURF'] = h*SURFc
# Biology
BIOc = ds['BLGBIOC']
budget['BIO'] = h*BIOc
# Correction and Forcing
CORRc = maskZ*ds['WTRAC01']/ds['drC']
budget['CORR'] = eta*CORRc
FORCc = ds['ForcTr01']
budget['FORC'] = eta*FORCc
# Pressure solver correction
epsilon = ds['oceFWflx']/rho0 + ds['WVEL'].isel(Zl=0) - xgrid.diff(ds['ETAN_snaps'],'T')/dt
budget['EPS'] = maskZ*epsilon*ds['TRAC01']

# Signs in closed budget
signs = {'TEND':-1,'ADV':1,'DIFF':1,'SURF':1,'BIO':1,'CORR':-1,'FORC':1,'EPS':-1}

# Residual
budget['RES'] = (signs['TEND']*budget['TEND']
                 + signs['ADV']*budget['ADV'] + signs['DIFF']*budget['DIFF'] 
                 + signs['SURF']*budget['SURF'] + signs['BIO']*budget['BIO'] 
                 + signs['CORR']*budget['CORR'] + signs['FORC']*budget['FORC']
                 + signs['EPS']*budget['EPS'])

# Transpose variables to be the same orientation
budget = budget.transpose('time', 'Z', 'YC', 'XC')
budget = budget.chunk({'time':1,'Z':52,'YC':588,'XC':int(2160/4)})

# Save
savedir = '/local/projects/bSOSE_carbon/budget-DIC/netcdf/'
for time in budget['time']:
    timepd = pd.to_datetime(time.values)
    outfile = 'bsose_i133_2013to2018_5day_'+str(timepd.date())+'_budget-DIC'
    path = savedir+outfile+'.nc'
    if os.path.isfile(path):
        if os.stat(path).st_size==5029554312:
            print('Already saved : '+outfile)
        else:
            print('Deleting partial file : '+outfile)
            os.system("rm ' + path")
            print('Saving : '+outfile)
            select = {'time':time}
            dsnow = budget.sel(select).expand_dims(dim='time')
            with ProgressBar():
                dsnow.to_netcdf(savedir+outfile+'.nc')
            dsnow.close()
    else:
        print('Saving : '+outfile)
        select = {'time':time}
        dsnow = budget.sel(select).expand_dims(dim='time')
        with ProgressBar():
            dsnow.to_netcdf(savedir+outfile+'.nc')
        dsnow.close()