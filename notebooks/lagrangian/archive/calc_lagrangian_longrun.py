# Master file for calculating trajectories in BSOSE

import xarray as xr
import pandas as pd
from xgcm import Grid
import numpy as np
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar
import os
from datetime import timedelta

from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D, ErrorCode, Variable

# Parameters
# ----------
N = 50

dt = 360 # minutes
outputdt = 30 # days
runtime = 30*365 # days

locinit = "Drake" # Drake, Ross
timedir = "back" # forw, back

# Returns
# -------
fileout = ("output"
           +".locinit_"+locinit
           +".timedir_"+timedir
           +".ntime_"+str(runtime)+".nc")

print("Saving trajectories to : "+fileout)

############

# SAMPLING KERNELS for T & S
# --------------------------
class TSMLDParticle(JITParticle):
    T = Variable('T', dtype=np.float32)
    S = Variable('S', dtype=np.float32)
    MLD = Variable('MLD', dtype=np.float32)

def SampleTSMLD(particle, fieldset, time):
    particle.T = fieldset.T[time,
                            particle.depth,
                            particle.lat,
                            particle.lon]
    particle.S = fieldset.S[time,
                            particle.depth,
                            particle.lat,
                            particle.lon]
    particle.MLD = fieldset.MLD[time,
                                particle.depth,
                                particle.lat,
                                particle.lon]

def mixedlayerdelete(particle, fieldset, time):
    # Absolute depth
    if -1*particle.depth <= particle.MLD:
        particle.delete()
        
def particledelete(particle, fieldset, time):
    particle.delete()
        
# TIME-LOOPING
# ------------
# Check if time-loop required based on runtime
if runtime>6*365:
    time_periodic = timedelta(runtime)
else:
    time_periodic=False

# FIELD SET
# ---------
rootdir = "/local/data/bSOSE/iter133NEW/5day/"
filenames = {'U':rootdir+'bsose_i133_2013to2018_5day_Uvel.nc',
             'V':rootdir+'bsose_i133_2013to2018_5day_Vvel.nc',
             'W':rootdir+'bsose_i133_2013to2018_5day_Wvel.nc',
             'T':rootdir+'bsose_i133_2013to2018_5day_Theta.nc',
             'S':rootdir+'bsose_i133_2013to2018_5day_Salt.nc',
             'MLD':rootdir+'bsose_i133_2013to2018_5day_MLD.nc'}
variables = {'U': 'UVEL',
             'V': 'VVEL',
             'W': 'WVEL',
             'T': 'THETA',
             'S': 'SALT',
             'MLD': 'BLGMLD'}
dimensions = {'U':{'lon':'XG','lat':'YC','depth':'Z','time':'time'},
             'V':{'lon':'XC','lat':'YG','depth':'Z','time':'time'},
             'W':{'lon':'XC','lat':'YC','depth':'Zl','time':'time'},
             'T':{'lon':'XC','lat':'YC','depth':'Z','time':'time'},
             'S':{'lon':'XC','lat':'YC','depth':'Z','time':'time'},
             'MLD':{'lon':'XC','lat':'YC','time':'time'}}
fs = FieldSet.from_netcdf(filenames,variables,dimensions,
                          deferred_load=True,gridindexingtype='mitgcm',
                          time_periodic=time_periodic)

# Add conditions for periodic boundary
fs.add_constant('halo_west', fs.U.grid.lon[0])
fs.add_constant('halo_east', fs.U.grid.lon[-1])
fs.add_periodic_halo(zonal=True)

def periodicBC(particle, fieldset, time):
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west

# PARTICLE SET
# ------------

# Get runtime details based on parameters
if locinit=="Drake":
    lats1D = np.linspace(-61.5,-55.5,N)
    lons1D = 291
    depths1D = -1*np.linspace(200,2200,N)
elif locinit=="Ross":
    lats1D = np.linspace(-75.5,-70,N)
    lons1D = 210
    depths1D= -1*np.linspace(100,3000,N)

elif locinit=="Drake-upstream":
    lats1D = np.linspace(-80,-55,N)
    lons1D = np.arange(231,291,10)
    depths1D = -1*np.linspace(200,3000,N)
    
[lons,lats,depths] = np.meshgrid(lons1D,lats1D,depths1D)

# Get start time based on direction
if timedir == "back":
    times=fs.U.grid.time[-1]
    dtsign = -1
elif timedir == "forw":
    times=fs.U.grid.time[0]
    dtsign = 1
    
pset = ParticleSet(fieldset=fs, pclass=TSMLDParticle,
                  lon=lons,lat=lats,depth=depths,time=times)

# Output save
if os.path.isfile(fileout):
    os.remove(fileout)
output_file = pset.ParticleFile(name=fileout, outputdt=timedelta(days=outputdt))

# RUN
# ---
kernel = AdvectionRK4_3D + pset.Kernel(SampleTSMLD) + pset.Kernel(mixedlayerdelete) + pset.Kernel(periodicBC)
pset.execute(kernel,
             runtime=timedelta(days=runtime),
             dt=dtsign*timedelta(minutes=dt),
             output_file=output_file,
             recovery={ErrorCode.ErrorOutOfBounds: particledelete})
output_file.export()
