import xarray as xr
from xgcm import Grid

def load_bsose():
    
    rootdir = '/local/data/bSOSE/'
    niter = 'iter133NEW'
    freq = '5day'
    chunks = 1

    ### CARBON CONCENTRATION BUDGET ###
    # Load all carbon tendencies and fluxes
    filenames = 'bsose_i133_2013to2018_5day_*C.nc'
    ds = xr.open_mfdataset(rootdir+niter+'/'+freq+'/'+filenames,chunks={'time':chunks})
    # Load surface carbon flux
    filenames = 'bsose_i133_2013to2018_5day_surfCO2flx.nc'
    ds_surf = xr.open_mfdataset(rootdir+niter+'/'+freq+'/'+filenames,chunks={'time':chunks})
    # Load carbon snapshots
    filenames = 'bsose_i133_2013to2018_5daySnapShots_DIC.nc'
    ds_csnaps = xr.open_dataset(rootdir+niter+'/'+freq+'/'+filenames,chunks={'time':chunks})
    ds_csnaps = ds_csnaps.rename({'time':'time_snaps','TRAC01':'TRAC01_snaps'}).drop('iter')
    # Specify shift of time axis for snapshots
    ds_csnaps['time_snaps'].attrs['c_grid_axis_shift']=-0.5

    ### VOLUME BUDGET ###
    # Load velocity data
    filenames = 'bsose_i133_2013to2018_5day_*vel.nc'
    ds_vel = xr.open_mfdataset(rootdir+niter+'/'+freq+'/'+filenames,chunks={'time':chunks})
    # Load SSH
    filenames = 'bsose_i133_2013to2018_5day_SSH.nc'
    ds_ssh = xr.open_mfdataset(rootdir+niter+'/'+freq+'/'+filenames,chunks={'time':chunks})
    # Load SSH snapshots
    filenames = 'bsose_i133_2013to2018_5daySnaps_SSH.nc'
    ds_sshsnaps = xr.open_mfdataset(rootdir+niter+'/'+freq+'/'+filenames,chunks={'time':chunks})
    ds_sshsnaps = ds_sshsnaps.rename({'time':'time_snaps','ETAN':'ETAN_snaps'}).drop('iter')
    # Specify shift of time axis
    ds_sshsnaps['time_snaps'].attrs['c_grid_axis_shift']=-0.5
    # Load FW flux
    freq = '1day'
    filenames = 'bsose_i133_2013to2018_1dy_oceFWflx.nc'
    ds_fw = xr.open_mfdataset(rootdir+niter+'/'+freq+'/'+filenames,chunks={'time':chunks})
    ds_fw = ds_fw.coarsen(time=5,boundary='trim',keep_attrs=True).mean().assign_coords({'time':ds_vel['time']})

    # Merge to full dataset
    ds = xr.merge([ds,ds_csnaps,ds_vel,ds_surf,ds_ssh,ds_sshsnaps,ds_fw])

    # Define vertical metrics as negative, to account for descending coordinate
    ds['drW'] = ds.hFacW * ds.drF #vertical cell size at u point
    ds['drS'] = ds.hFacS * ds.drF #vertical cell size at v point
    ds['drC'] = ds.hFacC * ds.drF #vertical cell size at tracer point
    # Volume
    ds['vC'] = ds['drC']*ds['rA']
    # Define cell side areas
    ds['rAW'] = ds['dyG']*ds['drW']
    ds['rAS'] = ds['dxG']*ds['drS']

    metrics = {
        ('X',): ['dxC', 'dxG'], # X distances
        ('Y',): ['dyC', 'dyG'], # Y distances
        ('Z',): ['drW', 'drS', 'drC'], # Z distances
        ('X', 'Y'): ['rA', 'rAs', 'rAw'] # Areas
    }

    xgrid = Grid(ds,periodic=['X'],metrics=metrics)
    
    return ds, xgrid

def create_xgcm_grid(ds):
    # Define vertical metrics as negative, to account for descending coordinate
    ds['drW'] = ds.hFacW * ds.drF #vertical cell size at u point
    ds['drS'] = ds.hFacS * ds.drF #vertical cell size at v point
    ds['drC'] = ds.hFacC * ds.drF #vertical cell size at tracer point
    # Volume
    ds['vC'] = ds['drC']*ds['rA']
    # Define cell side areas
    ds['rAW'] = ds['dyG']*ds['drW']
    ds['rAS'] = ds['dxG']*ds['drS']

    metrics = {
        ('X',): ['dxC', 'dxG'], # X distances
        ('Y',): ['dyC', 'dyG'], # Y distances
        ('Z',): ['drW', 'drS', 'drC'], # Z distances
        ('X', 'Y'): ['rA', 'rAs', 'rAw'] # Areas
    }

    xgrid = Grid(ds,periodic=['X'],metrics=metrics)
    return ds, xgrid
    