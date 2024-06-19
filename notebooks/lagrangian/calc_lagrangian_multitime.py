#!/home/aos/graemem/miniconda3/envs/parcels/bin/python
# Master file for calculating trajectories in BSOSE

import xarray as xr
import pandas as pd
from xgcm import Grid
import numpy as np
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar
import os
from datetime import timedelta

from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D, StatusCode, Variable

# Parameters
# ----------
Nh = 10
Nz = 20

dt = 360 # minutes
outputdt = 10 # days
runtime = 6*365 #1*365 # days

locinit = "Shelf" # Drake, Drake-upstream, Ross, coral, coral-approx, RossShelf.?, Shelf.?
timedir = "back" # forw, back
multitimes = False # initialize at multiple different times
Nt = 12 # Number of initialized times
Dt = 6 # 5-day multiples; spacing between initialization time

# For locinit==shelf, locations are derived in calc_init_shelf.ipynb and loaded from there

# Returns
# -------
fileout = ("output"
           +".locinit_"+locinit
           +".Nh_"+str(Nh)
           +".Nz_"+str(Nz)
           +".timedir_"+timedir
           +".ntime_"+str(runtime)
           +".dt_"+str(dt))
if multitimes:
    fileout+=".itimes_N"+str(Nt)+"_D"+str(Dt)
fileout+=".zarr"

print("Saving trajectories to : "+fileout)

############

# KERNELS
# --------------------------
# SAMPLING T, S, MLD
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

# Deleting particles that enter the mixed layer
def mixedlayerdelete(particle, fieldset, time):
    # Absolute depth
    if -1*particle.depth <= particle.MLD:
        particle.delete()
        
def particledelete(particle, fieldset, time):
    particle.delete()

def DeleteErrorParticle(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()
        
# TIME-LOOPING
# ------------
# Check if time-loop required based on runtime
if multitimes:
    maxt = runtime+Nt*Dt*5
else:
    maxt = runtime
    
if maxt>=6*365:
    time_periodic = timedelta(runtime)
else:
    time_periodic=False

# FIELD SET
# ---------
rootdir = "/gws/nopw/j04/co2clim/datasets/bSOSE/ITER133/"
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
        particle_dlon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle_dlon -= fieldset.halo_east - fieldset.halo_west

# PARTICLE SET
# ------------

# Get runtime details based on parameters
if locinit=="Drake":
    lats1D = np.linspace(-61.5,-55.5,Nh)
    lons1D = 291
    depths1D = -1*np.linspace(200,2200,Nh)
elif locinit=="Ross":
    lats1D = np.linspace(-75.5,-70,Nh)
    lons1D = 210
    depths1D= -1*np.linspace(100,3000,Nz)
elif locinit=="Drake-upstream":
    lats1D = np.linspace(-80,-55,Nh)
    lons1D = np.arange(201,291,10)
    depths1D = -1*np.linspace(200,3000,Nz)
else:
    lats1D = np.empty(0)
    lons1D = np.empty(0)
    depths1D = np.empty(0)
    
[lons,lats,depths] = np.meshgrid(lons1D,lats1D,depths1D)

# To initialize at approximate coral locations
if locinit=="coral":
    df = pd.read_excel('Coral_Locations.xlsx',engine='openpyxl')
    lats = np.empty(0)
    lons = np.empty(0)
    depths = np.empty(0)
    for name in df['names'].unique():
        lon = 360-df[df['names']==name]['Longitude'].iloc[0]
        lat = -1*df[df['names']==name]['Latitude'].iloc[0]
        depthmin = -1*df[df['names']==name]['Depth'].min()
        depthmax = -1*df[df['names']==name]['Depth'].max()

        latV = [lat-0.5,lat+0.5]
        depthV = [np.round(depthmax,-1)-10,np.round(depthmin,-1)+10]

        lons1D = lon
        lats1D = np.linspace(latV[0],latV[1],Nh)
        depths1D = np.linspace(depthV[0],depthV[1],Nz)

        [lons1,lats1,depths1] = np.meshgrid(lons1D,lats1D,depths1D)

        lons = np.concatenate((lons,lons1.flatten()))
        lats = np.concatenate((lats,lats1.flatten()))
        depths = np.concatenate((depths,depths1.flatten()))
elif locinit=="coral-approx":
    locs = {'pink':{'lat':[-60.5,-59.5],'depth':[-2000,-400]},
           'green':{'lat':[-57.5,-56.5],'depth':[-2000,-400]}}
    lats = np.empty(0)
    lons = np.empty(0)
    depths = np.empty(0)
    for loc in locs.keys():
        lon = 360-66
        latV = locs[loc]['lat']
        depthV = locs[loc]['depth']

        lons1D =  lon
        lats1D = np.linspace(latV[0],latV[1],Nh)
        depths1D = np.arange(depthV[0],depthV[1]+10,10)

        [lons1,lats1,depths1] = np.meshgrid(lons1D,lats1D,depths1D)

        lons = np.concatenate((lons,lons1.flatten()))
        lats = np.concatenate((lats,lats1.flatten()))
        depths = np.concatenate((depths,depths1.flatten()))

# To initialise on Antarctic shelf
if (locinit=="RossShelf") or (locinit=="Shelf"):
    if locinit=="RossShelf":
        path = "../../data/RossShelf.level1000.sep500.maxd5000.txt"
    elif locinit=="Shelf":
        path = "../../data/Shelf.level-1000.sep-500.txt"
    init = np.loadtxt(path)
    # hack to correct for different way that initialisation locations are saved
    if locinit=="Shelf":
        init=np.transpose(init)
    initlon = init[0,:]
    initlat = init[1,:]
    dh = 0.1
    
    lats = np.empty(0)
    lons = np.empty(0)
    depths = np.empty(0)

    for i in range(len(initlon)):
        lon = initlon[i]
        lat = initlat[i]

        depthmin = 0
        depthmax = 1000

        lonV = [lon-dh/2,lon+dh/2]
        latV = [lat-dh/2,lat+dh/2]
        depthV = [depthmin,depthmax]

        lons1D =  np.linspace(lonV[0],lonV[1],Nh)
        lats1D = np.linspace(latV[0],latV[1],Nh)
        depths1D = -1*np.linspace(depthV[0],depthV[1],Nz)

        [lons1,lats1,depths1] = np.meshgrid(lons1D,lats1D,depths1D)

        lons = np.concatenate((lons,lons1.flatten()))
        lats = np.concatenate((lats,lats1.flatten()))
        depths = np.concatenate((depths,depths1.flatten()))
        
# Get start time based on direction
if timedir == "back":
    itimes = -1
    dtsign = -1
elif timedir == "forw":
    itimes = 0
    dtsign = 1
    
if multitimes:
    Dt = dtsign*Dt
    ftime = itimes
    ltime = ftime+Dt*Nt
    itimes = np.arange(ftime,ltime,Dt)
    
times=fs.U.grid.time[itimes]

if multitimes:
    npart = len(lons.flatten())
    times = np.repeat(times,npart)
    
    lons = np.tile(lons.flatten(),Nt)
    lats = np.tile(lats.flatten(),Nt)
    depths = np.tile(depths.flatten(),Nt)

print("ntimes: "+str(len(times.flatten())))
print("nlons: "+str(len(lons.flatten())))
print("nlats: "+str(len(lats.flatten())))
print("ndepths: "+str(len(depths.flatten())))

pset = ParticleSet(fieldset=fs, pclass=TSMLDParticle,
                  lon=lons,lat=lats,depth=depths,time=times)

# Output save
if os.path.isfile(fileout):
    os.remove(fileout)
output_file = pset.ParticleFile(name=fileout, outputdt=timedelta(days=outputdt))

# RUN
# ---
kernel = AdvectionRK4_3D + pset.Kernel(SampleTSMLD) + pset.Kernel(mixedlayerdelete) +  pset.Kernel(periodicBC) + pset.Kernel(DeleteErrorParticle)
pset.execute(kernel,
             runtime=timedelta(days=runtime),
             dt=dtsign*timedelta(minutes=dt),
             output_file=output_file)
