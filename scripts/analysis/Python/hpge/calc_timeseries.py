#!/usr/bin/env python

#     |------------------------------------------------------------|
#     | This module creates a 2D field of maximum spurious current |
#     | in the vertical and in time after an HPGE test.            |
#     | The resulting file can be used then to optimise the rmax   |
#     | of Multi-Envelope vertical grids.                          |
#     |                                                            |
#     | Author: Diego Bruciaferri                                  |
#     | Date and place: 07-09-2021, Met Office, UK                 |
#     |------------------------------------------------------------|


import os
from os.path import join, isfile, basename, splitext
import glob
import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar
from utils import compute_masks

# ==============================================================================
# Input files
# ==============================================================================

# Folder path containing HPGE spurious currents velocity files 
MAINdir = '/scratch/dbruciaf/SE-NEMO/hpge/'
HPGElst = 'u-cg602_hpge_se-nemo_r018-010-010_glo-r018-010_ant_opt_v3_3months_traldfoff'
DOMCFG = '/data/users/dbruciaf/SE-NEMO/se-orca025/se-nemo-domain_cfg/domain_cfg_MEs.nc'
loc_msk_glo = '/data/users/dbruciaf/SE-NEMO/se-orca025/MEs_450-800_3200/bathymetry.loc_area.dep3200_polglo_sig3_itr1.nc'
loc_msk_ant = '/data/users/dbruciaf/SE-NEMO/se-orca025/MEs_450-800_3200/bathymetry.loc_area.dep3200_polant_sig3_itr1.nc'

label = 'MEs'

# Name of the zonal and meridional velocity variables
Uvar = 'uo'
Vvar = 'vo'
# Name of the variable to chunk with dask and size of chunks
chunk_var = 'time_counter'
chunk_size = 1

# ==============================================================================
# loc msk
ds_loc_glo  = xr.open_dataset(loc_msk_glo).squeeze()
msk_glo = ds_loc_glo['s2z_msk']
msk_glo = msk_glo.where(msk_glo==0,1)
ds_loc_ant  = xr.open_dataset(loc_msk_ant).squeeze()
msk_ant = ds_loc_ant['s2z_msk']
msk_ant = msk_ant.where(msk_ant==0,1)

# Loading domain geometry
ds_dom  = xr.open_dataset(DOMCFG).squeeze()

# Computing land-sea masks
ds_dom = compute_masks(ds_dom, merge=True)

e3t = ds_dom["e3t_0"].squeeze()
e2t = ds_dom["e2t"].squeeze()
e1t = ds_dom["e1t"].squeeze()

e1t = e1t.where(ds_dom.tmask==1)
e2t = e2t.where(ds_dom.tmask==1)
e3t = e3t.where(ds_dom.tmask==1)

e1t = e1t.where((msk_glo==1) & (msk_ant==1))
e2t = e2t.where((msk_glo==1) & (msk_ant==1))
e3t = e3t.where((msk_glo==1) & (msk_ant==1))

cel_vol = e1t * e2t * e3t
dom_vol = cel_vol.sum(skipna=True)

HPGEdir = MAINdir + HPGElst

Ufiles = sorted(glob.glob(HPGEdir+'/*grid_U*.nc'))
Vfiles = sorted(glob.glob(HPGEdir+'/*grid_V*.nc'))

v_max     = []
v_99p_tot = []
v_99p_loc = []
v_avg     = []

for F in range(len(Ufiles)):

    print(Ufiles[F])

    ds_U = xr.open_dataset(Ufiles[F], chunks={chunk_var:chunk_size})
    ds_V = xr.open_dataset(Vfiles[F], chunks={chunk_var:chunk_size})
    U4   = ds_U[Uvar]
    V4   = ds_V[Vvar]

    # rename some dimensions
    U4 = U4.rename({U4.dims[0]: 't', U4.dims[1]: 'z'})
    V4 = V4.rename({V4.dims[0]: 't', V4.dims[1]: 'z'})

    # interpolating from U,V to T
    U = U4.rolling({'x':2}).mean().fillna(0.)
    V = V4.rolling({'y':2}).mean().fillna(0.) 
    
    vel_t = np.sqrt(np.power(U,2) + np.power(V,2)) 
    vel_l = vel_t.where(ds_dom.tmask==1)
    vel_l = vel_l.where((msk_glo==1) & (msk_ant==1))

    v_max.extend(vel_t.max(dim=('z','y','x')).values.tolist())
    v_99p_tot.extend(vel_t.quantile(0.99,dim=('z','y','x')).values.tolist())
    v_99p_loc.extend(vel_l.quantile(0.99,dim=('z','y','x')).values.tolist())
    v_avg.extend(((cel_vol*vel_l).sum(dim=["x","y","z"], skipna=True) / dom_vol).values.tolist())

# Saving 

ds = xr.Dataset()
ds["max_u"] = xr.DataArray(np.asarray(v_max), dims=('t'))
ds["u_99p_tot"] = xr.DataArray(np.asarray(v_99p_tot), dims=('t'))
ds["u_99p_loc"] = xr.DataArray(np.asarray(v_99p_loc), dims=('t'))
ds["avg_u"] = xr.DataArray(np.asarray(v_avg), dims=('t'))

# -------------------------------------------------------------------------------------   
# Writing the max_hpge file

print('WRITING the maximum_hpge.nc FILE')

out_file = "hpge_timeseries.nc"
delayed_obj = ds.to_netcdf(join(HPGEdir,out_file), compute=False)

with ProgressBar():
     results = delayed_obj.compute()

