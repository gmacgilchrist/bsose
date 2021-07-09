import os,sys
import numpy as np
import netCDF4 as nc
import matplotlib as mpl

import matplotlib.pyplot as plt
from matplotlib import patches as patches
from scipy.io import loadmat,savemat
from scipy import interpolate

from datetime import datetime

# run clustering algorithms for T/S, using PCA as dimension reduction
from sklearn.preprocessing import Imputer	
from sklearn import preprocessing#, mixture
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.mixture import GaussianMixture as GMM
from sklearn.model_selection import train_test_split

# import colormaps
from palettable.colorbrewer.diverging import PRGn_10 as prg

#~~~~~~~~~~~~~~~~~~~~~~~~
# time interval
tsnap     = 1 
tmax      = 12*5 # only 2008 
dt        = 2628900. # monthly = 86400 * 30.427083333333332
nt        = tmax-tsnap+1
#TimeStep  = 24:24:43848;

# folder to save the files in
#dir_path = '/data/irosso/data/BSOSE/DIC/3D/DAILY/200pts_6/'

# files to read 
maindir    = '/data/SOSE/SOSE/SO6/ITER122'
bgtdir     = '/data/SOSE/SOSE/SO6/ITER122/budgets'

# MITgcm grid
grid_dir   = '/data/soccom/GRID_6/grid.nc'
data       = nc.Dataset(grid_dir)
XC         = data.variables['XC'][:]
YC         = data.variables['YC'][:]
RC         = data.variables['RC'][:]
RAC        = data.variables['RAC'][:]
DXG        = data.variables['DXG'][:]
DYG        = data.variables['DYG'][:]
DRF        = data.variables['DRF'][:]
hFacW      = data.variables['hFacW'][:]
hFacS      = data.variables['hFacS'][:]
hFacC      = data.variables['hFacC'][:]

# size grid
len_x      = len(XC[0,:])
len_y      = len(YC[:,0])
len_z      = len(RC)

# extract only one portion
xel        = np.arange(0,100,1)#80,1)
yel        = np.arange(0,len(YC[:,0])-1,1)     
zel        = np.arange(0,51,1)

nx         = len(xel)
ny         = len(yel)
nz         = len(zel)

# cell volume, face areas (for flux calculations)
volume     = np.zeros((nz,ny,nx),'>f4')
AREAWEST   = np.zeros((nz,ny,nx),'>f4')
AREASOUTH  = np.zeros((nz,ny,nx),'>f4')
AREACELL   = np.zeros((nz,ny,nx),'>f4')
for k in range(nz):
    volume[k,...]    = (hFacC[k,yel[0]:yel[-1]+1,xel]*RAC[yel[0]:yel[-1]+1,xel].transpose()*DRF[k]).transpose()
    AREACELL[k,...]  = RAC[yel[0]:yel[-1]+1,xel]
    AREAWEST[k,...]  = DYG[yel,xel[0]:xel[-1]+1]*DRF[k]*hFacW[k,yel,xel[0]:xel[-1]+1]
    AREASOUTH[k,...] = DXG[yel[0]:yel[-1]+1,xel]*DRF[k]*hFacS[k,yel[0]:yel[-1]+1,xel].transpose()

DZ         = np.nan*np.ones((nz,ny,nx),'>f4')
DRF1       = np.nan*np.ones((nz,ny,nx),'>f4')
DRF        = DRF[:nz]
for kk in range(nz):
	DRF1[kk,...] = DRF[kk]
	DZ[kk,...]   = DRF[kk]*hFacC[k,yel[0]:yel[-1]+1,xel].transpose()

# read snapshots
file       = os.path.join(bgtdir,'bsose_i122_2013to2017_MonthlySnapShots_Fe.nc') # in mol Fe
data       = nc.Dataset(file)
Tr         = data.variables['TRAC06'][:]

# tendency of tracer 
# (change in tracer content / storage)
dTrdt      = ( Tr[1:,zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1] -  Tr[:-1,zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1])/dt
# save the file 
#file_name = os.path.join(dir_path, 'dFedt_BL2_3D.mat')
#savemat(file_name,'dTrdt', '-v7.3').....

# prepare the arrays for the bgt
# Advection u.(nabla*Tr)
ADV_x      = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
ADV_y      = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
ADV_z      = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
# Divergence (nabla.u) * tr
DIV_x      = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
DIV_y      = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
DIV_z      = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
# Diffusion
DIFF_x     = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
DIFF_y     = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
DIFF_z     = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
# Biology
BIO        = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
# Sediment
SED        = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
# Aeolian Deposition
DEPO       = np.nan*np.zeros((nt,nz,ny,nx),'>f4')
# Correction due to vertical advection (WTRAC06)
CORR       = np.zeros((nt,nz,ny,nx),'>f4')

# load the other terms, timestep by timestep
for tt in range(nt):
	file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthlyADVx_FE.nc')
	data       = nc.Dataset(file)
	FLUXx      = data.variables['ADVxTr06'][tt,:,:,:]

	file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthlyADVy_FE.nc')
	data       = nc.Dataset(file)
	FLUXy      = data.variables['ADVyTr06'][tt,:,:,:]

	file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthlyADVr_FE.nc')
	data       = nc.Dataset(file)
	FLUXz      = data.variables['ADVrTr06'][tt,:,:,:]

	# tendency due to biology
	# NOTE: this includes flux from sediments
	# - uptake + remin + sed (mol/m3/s)
	file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthly_BLGBIOFE.nc')
	data       = nc.Dataset(file)
	BIOf       = data.variables['BLGBIOFE'][tt,:,:,:]

	file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthly_BLGFESED.nc')
	data       = nc.Dataset(file)
	SEDf       = data.variables['BLGFESED'][tt,:,:,:]
	BIOf       = BIOf-SEDf

	file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthlyDFxE_FE.nc')
	data       = nc.Dataset(file)
	DIFFx      = data.variables['DFxETr06'][tt,:,:,:]

	file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthlyDFyE_FE.nc')
	data       = nc.Dataset(file)
	DIFFy      = data.variables['DFyETr06'][tt,:,:,:] 

	file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthlyDFrI_FE.nc')
	data       = nc.Dataset(file)
	DIFFz      = data.variables['DFrITr06'][tt,:,:,:]

	file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthlyForcFE.nc') # is this the surface flux (aeolian deposition)???
	data       = nc.Dataset(file)
	FORC       = data.variables['ForcTr06'][tt,:,:,:]

	#file       = os.path.join(bgtdir,'bsose_i122_2013to2017_monthly_WFE.nc')
	#data       = nc.Dataset(file)
	#fieldCorr        = data.variables['WTRAC06'][tt,:,:,:]

	file       = os.path.join(maindir,'bsose_i122_2013to2017_MonthlySnapShots_UVEL.nc')
	data       = nc.Dataset(file)
	fieldU     = data.variables['UVEL'][tt,:,:,:]

	file       = os.path.join(maindir,'bsose_i122_2013to2017_MonthlySnapShots_VVEL.nc')
	data       = nc.Dataset(file)
	fieldV     = data.variables['VVEL'][tt,:,:,:]

	file       = os.path.join(maindir,'bsose_i122_2013to2017_MonthlySnapShots_WVEL.nc')
	data       = nc.Dataset(file)
	fieldW     = data.variables['WVEL'][tt,:,:,:]
	
	file	   = os.path.join(maindir,'bsose_i122_2013to2017_monthly_WFE.nc')
	data       = nc.Dataset(file)
	fieldCorr  = data.variables['WTRAC06'][tt,:,:,:]

	# Advection (div(bar(u)*tr))
	ADV_x[tt,...] = np.diff(FLUXx[zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+2],axis=2)/volume
	ADV_y[tt,...] = np.diff(FLUXy[zel[0]:zel[-1]+1,yel[0]:yel[-1]+2,xel[0]:xel[-1]+1],axis=1)/volume
	ADV_z[tt,...] = -np.diff(FLUXz[zel[0]:zel[-1]+2,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1],axis=0)/volume

	# terms for the equation 
	ADVn        = -(ADV_x[tt,...]+ADV_y[tt,...]+ADV_z[tt,...])

	# save the terms
	#file_name = os.path.join(dir_path, 'ADV_Fe_BL2_3D.mat')
	#savemat(file_name,'ADVn','ADV_x','ADV_y','ADV_z','-v7.3')

	# Divergence (nabla.u) * tr
    U     = fieldU[zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1]*AREAWEST
    V     = fieldV[zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1]*AREASOUTH
    W     = fieldW[zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1]*AREACELL
    W[0,...] = 0

    # gradients:
    DIV_z[tt,:-1,...]  = ((W[1:,:,:]-W[:-1,:,:])/volume[:-1,...])*Tr[tt,zel[0]:zel[-1],yel[0]:yel[-1]+1,xel[0]:xel[-1]+1]
    DIV_y[tt,:,:-1,:]  = ((V[:,1:,:]-V[:,:-1,:])/volume[:,:-1,:])*Tr[tt,zel[0]:zel[-1]+1,yel[0]:yel[-1],xel[0]:xel[-1]+1]
    DIV_x[tt,...,:-1]  = ((U[:,:,:-1]-U[:,:,1:])/volume[...,:-1])*Tr[tt,zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]]
      
    # compute DIC*DIV_vec(U)
    DIV_t  = DIV_x[tt,...]+DIV_y[tt,...]+DIV_z[tt,...]

	# save the terms
	#file_name = os.path.join(dir_path, 'DIV_Fe_BL2_3D.mat')
	#save([file_name],'DIV_t','DIV_x','DIV_y','DIV_z','-v7.3');

	# Diffusion
	DIFF_x[tt,...]  = np.diff(DIFFx[zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+2],axis=2)/volume
	DIFF_y[tt,...]  = np.diff(DIFFy[zel[0]:zel[-1]+1,yel[0]:yel[-1]+2,xel[0]:xel[-1]+1],axis=1)/volume
	DIFF_z[tt,...]  = -np.diff(DIFFz[zel[0]:zel[-1]+2,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1],axis=0)/volume

	DIFF_horiz  = DIFF_x[tt,...] + DIFF_y[tt,...]
	DIFF        = DIFF_horiz + DIFF_z[tt,...] 
	DIFFn       = -DIFF

	# save the terms
	#file_name = os.path.join(dir_path, 'DIFF_Fe_BL2_3D.mat')
	#save([file_name],'DIFFn', 'DIFF_horiz','DIFF_z','-v7.3')

	# Biological term
	BIO[tt,...] = BIOf[zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1]
	# save the terms
	#file_name  =  [dir_path 'BIO_BL2_3D.mat']
	#save([file_name],'BIO','-v7.3')

	# Sediment term
	SED[tt,...] = SEDf[zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1]
	# save the terms
	#file_name  =  [dir_path 'BIO_BL2_3D.mat']
	#save([file_name],'BIO','-v7.3')

	# Aeolian deposition
	DEPO[tt,...] = FORC[zel[0]:zel[-1]+1,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1]

	# correction to vertical advection at z=0
	CORR[tt,0,...] = fieldCorr[0,yel[0]:yel[-1]+1,xel[0]:xel[-1]+1]/DZ[0,...]

	# remove correction from advection:
	ADVcorr = ADVn[tt,...] - CORR[tt,...]
	
	RES = -dTrdt[tt,...] + ADVcorr + DIFFn + BIO[tt,...] + SED[tt,...]  + DEPO[tt,...] + DIV_t 
	
	plt.plot(RES[:,400,50],RC[:-1],'k--',label='res')
	plt.plot(dTrdt[tt,:,400,50],RC[:-1],'k',label='tendency')
	plt.plot(ADVcorr[:,400,50],RC[:-1],'b',label='-adv')
	plt.plot(DIFFn[:,400,50],RC[:-1],'red',label='-diff')
	plt.plot(BIO[tt,:,400,50],RC[:-1],'g',label='-bio')
	plt.plot(SED[tt,:,400,50],RC[:-1],'orange',label='sediment')
	plt.plot(DEPO[tt,:,400,50],RC[:-1],'brown',label='deposition')
	plt.legend(loc=1)
	plt.show()
	
	
