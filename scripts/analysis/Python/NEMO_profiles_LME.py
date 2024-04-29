#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Companion script to EN4_profiles_LME.py

Extracts profiles from LME regions that correspond to the profiles already extracted by
EN4_profiles_LME.py

Created on Tue Apr 25 13:41:16 2023

ssh sci6 -L 6000:localhost:22
pycharm jelt@loclahost:6000 key
env: senemo-profile-lme (over ssh)
COAsT branch:  feature/535_stratification_diag

@author: jelt
"""
import numpy as np
import matplotlib.pylab as plt
import scipy.io
import pandas as pd
import re
import sys
import time
import xarray as xr
import os
import socket
isliv = 'livljobs' in socket.gethostname()
if isliv:
 sys.path.insert(0,'/login/$USER/work/GitHub/COAsT/')
 #sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jelt/GitHub/COAsT/')

import coast


def extract_profiles_from_gridded(profile, gridded):
    """
    input: profile object containing indices on gridded object
            gridded object containing ...
    output: extracted_gridded_profiles, profile_matched - new synthetic profiles and input profiles, trimmed to match
    """
    extracted = coast.Profile(config=fn_profile_config)
    varlist = ["conservative_temperature", "absolute_salinity", "depth", "e3_0"]
    store2d = np.ones((len(varlist), profile.dataset.dims["id_dim"], gridded.dataset.dims["z_dim"] ))
    store1d = np.ones((3, profile.dataset.dims["id_dim"]))
    time_max = gridded.dataset.time.max().values
    time_min = gridded.dataset.time.min().values
    for count in range(profile.dataset.dims["id_dim"]):
        tI = nearest_time_index(gridded.dataset.time, profile.dataset.time.isel(id_dim=count))
        if (profile.dataset.time.isel(id_dim=count) <= time_max) & (profile.dataset.time.isel(id_dim=count) >= time_min):
            for ivar, var in enumerate(varlist):
                store2d[ivar, count, :] = gridded.dataset.isel(x_dim=int(profile.dataset.ind_i[count])).isel(y_dim=int(profile.dataset.ind_j[count])).isel(t_dim=tI)[var]
            store1d[0, count]  = gridded.dataset.time.isel(t_dim=tI)
            store1d[1, count] = gridded.dataset.isel(x_dim=int(profile.dataset.ind_i[count])).isel(
                y_dim=int(profile.dataset.ind_j[count]))['latitude']
            store1d[2, count] = gridded.dataset.isel(x_dim=int(profile.dataset.ind_i[count])).isel(
                y_dim=int(profile.dataset.ind_j[count]))['longitude']

        else:
            store2d[:, count, :] = np.nan
            store1d[:, count] = np.nan
            print('out of time range - skip')

    # Rebuild dataset as a Profile object
    da_time = xr.DataArray(data=store1d[0,:].flatten(), dims=["id_dim"]).astype("datetime64[ns]")
    da_lat = xr.DataArray(data=store1d[1,:].flatten(), dims=["id_dim"])
    da_lon = xr.DataArray(data=store1d[2,:].flatten(), dims=["id_dim"])

    da_CT  = xr.DataArray(data=store2d[0, :, :].squeeze(), dims=["id_dim", "z_dim"])
    da_AS  = xr.DataArray(data=store2d[1,:,:].squeeze(), dims=["id_dim", "z_dim"])
    da_dep = xr.DataArray(data=store2d[2,:,:].squeeze(), dims=["id_dim", "z_dim"])
    da_e3  = xr.DataArray(data=store2d[3,:,:].squeeze(), dims=["id_dim", "z_dim"])

    ds = xr.Dataset({"conservative_temperature": da_CT,
                     "absolute_salinity": da_AS,
                     "depth": da_dep,
                     "e3_0": da_e3,
                     "latitude": da_lat, "longitude": da_lon,
                     "time": da_time})

    keep_indices = ~np.isnan(ds.time).values
    ds = ds.isel(id_dim=keep_indices)
    profile_new = profile.isel(id_dim=keep_indices)

    #ds = ds.dropna(dim="id_dim", how="all")  # works but can not be resused on profile_new
    ds = ds.set_coords(["latitude", "longitude", "time"])

    extracted.dataset = ds
    return extracted, profile_new

def nearest_time_index(items, pivot):
    """ find the index for the closest timestamp """
    #return np.where(gridded.dataset.time == min(gridded.dataset.time, key=lambda x: abs(x - profile.dataset.time.isel(id_dim=count)) ))
    return np.where(items == min(items, key=lambda x: abs(x - pivot) ))[0][0].astype(int)

# load LME data for names
LME_Data=np.load('/home/users/jelt/GitHub/SE-NEMO/scripts/analysis/Data/LME_gridinfo_V4.npz')

ystart = 1978
ystop = 2019

month = 2
year = 1983
yr_str = str(year)
mm_str = str(month).zfill(2)

# Load the rebinned EN4 data by LME
LME_en4_in_dir = "/gws/nopw/j04/class_vol2/senemo/jelt/PROCESSED/EN4.2.1/1978-2019/"
in_path = LME_en4_in_dir
out_path = LME_en4_in_dir

#%% load model data
root = "../../../../"
# And by defining some file paths
#dn_files = root + "./example_files/"
#fn_nemo_dat = os.path.join(dn_files, "coast_example_nemo_data.nc")
#fn_nemo_dom = os.path.join(dn_files, "coast_example_nemo_domain.nc")
#fn_nemo_config = os.path.join(root, "./config/example_nemo_grid_t.json")

# Path to NEMO model
fn_nemo_data = "/gws/nopw/j04/class_vol2/senemo/FINAL_PHYSICS_RUNS/GS1p0_notide/output/SENEMO_1m_*grid_T*" + yr_str + mm_str + ".nc"
#fn_nemo_data = "/gws/nopw/j04/class_vol2/senemo/FINAL_PHYSICS_RUNS/GS1p0_notide/output/SENEMO_1m_*1983*grid_T*.nc"
#fn_nemo_data = "/gws/nopw/j04/class_vol2/senemo/FINAL_PHYSICS_RUNS/GS1p0_notide/output/SENEMO_1m*grid_T*.nc"
fn_nemo_grid = "/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/config/domain_cfg.nc"
fn_nemo_config_t = "/gws/nopw/j04/class_vol2/senemo/jelt/PROCESSED/tmp_nemo_grid_t.json"

fn_profile_config='/home/users/jelt/GitHub/COAsT/config/example_en4_profiles.json'
# Create gridded object:
nemo = coast.Gridded(fn_nemo_data, fn_nemo_grid, multiple=True, config=fn_nemo_config_t)

'''
#### Create a landmask array in Gridded
In this example we add a `landmask` variable to the `Gridded` dataset.
When this is present, the `obs_operator` will use this to interpolation to the
nearest *wet* point. If not present, it will just take the nearest grid point (not implemented).

We also rename the depth at initial time coordinate `depth_0` to `depth` as this is expected by Profile()
'''
nemo.dataset["landmask"] = nemo.dataset.bottom_level == 0
nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})  # profile methods will expect a `depth` coordinate


nLME = LME_Data["DOMNAM"].shape[0]
#for iLME in range(34,36): #66):
for iLME in range(66):

    LME_Name=LME_Data["DOMNAM"][iLME]

    start_time = time.perf_counter()
    print(LME_Name)

    fn_in = '{0}/{1}_{2}_{3}_EN4_PEA_SST_SSS_binned.nc'.format(in_path, LME_Name, ystart, ystop)

    profile = coast.Profile(config=fn_profile_config)
    profile.dataset = xr.open_dataset(fn_in) #, chunks={'id_dim': 10000})

    #profile = coast.Profile(config=fn_profile_config)
    #profile.dataset = pprofile.dataset.isel(id_dim=slice(0,4))
    #model_profiles = nemo.dataset.isel(x_dim=profile.dataset.ind_i).isel(y_dim=profile.dataset.ind_j)
    #profile.dataset = profile.dataset.sortby(profile.dataset.time)
    #profile.dataset = profile.dataset.expand_dims(dim={"z_dim": 1})
    #profile.dataset['depth'] = xr.DataArray(np.zeros((profile.dataset.dims["id_dim"], 1)), dims = ("id_dim", "z_dim"))

    # For the global model the wrapped seam can cause problems with the nearest neighbour calculation if the NEMO seam is within the LME
    #if ((nemo.dataset.longitude.isel(x_dim=0).mean() < profile.dataset.longitude.values.max())
    #    & (nemo.dataset.longitude.isel(x_dim=0).mean() > profile.dataset.longitude.values.min())):
    #    print(f"Need to roll the dataset {LME_Name}")
    #    nemo.dataset = nemo.dataset.roll(x_dim=nemo.dataset.dims["x_dim"]//2, roll_coords=True)

    # Use obs operator for horizontal remapping of Gridded onto Profile.
    #model_profiles = profile.obs_operator(nemo)
    model_profiles, profile = extract_profiles_from_gridded(profile, nemo)
    model_profiles.dataset['interp_dist'] = 0 * model_profiles.dataset.latitude
    model_profiles.dataset['interp_lag'] = model_profiles.dataset.time - model_profiles.dataset.time

    if(0):
        plt.scatter(model_profiles.dataset.longitude,
                    model_profiles.dataset.latitude,
                    c=model_profiles.dataset.interp_dist)
        plt.colorbar()
        plt.title(LME_Name)
        #plt.savefig(f'LOGS/{LME_Name}_1.png')
        plt.show()

        plt.plot( model_profiles.dataset.interp_dist, model_profiles.dataset.interp_lag.astype(np.timedelta64)/np.timedelta64(1, 'D'),'.')
        plt.xlabel('interp_dist (km)')
        plt.ylabel('inter_lag (days)')
        #plt.savefig(f'LOGS/{LME_Name}_2.png')
        plt.show()

        #plt.scatter(
        #    profile.dataset.time.diff("id_dim").astype(np.timedelta64) / np.timedelta64(1, 'D'),
        #    c=profile.dataset.interp_dist[0:-2]);
        #plt.colorbar();
        #plt.show()

        plt.plot(model_profiles.dataset.time, model_profiles.dataset.interp_dist,'.')
        plt.xlabel('time')
        plt.ylabel('interp_dist (km)')
        #plt.savefig(f'LOGS/{LME_Name}_3.png')
        plt.show()

    if(1):
        ### Discard profiles where the interpolation distance is too large
        too_far = 13.75  # 110/4/2 distance km
        keep_indices = model_profiles.dataset.interp_dist <= too_far
        model_profiles = model_profiles.isel(id_dim=keep_indices)
        profile = profile.isel(id_dim=keep_indices)
        if np.sum(~keep_indices.values) > 0:
            print(f"Dropped {np.sum(~keep_indices.values)} of {len(keep_indices)} profiles: too far in space")

        # Throw away profile where the interpolation time is larger than 48h
        keep_indices = np.abs(model_profiles.dataset.interp_lag) <= np.timedelta64(48,
                                                                                   'h')  ## SHOULD MOVE PARAMTER TO CONFIG
        model_profiles = model_profiles.isel(id_dim=keep_indices)
        profile = profile.isel(id_dim=keep_indices)
        if np.sum(~keep_indices.values) > 0:
            print(f"Dropped {np.sum(~keep_indices.values)} profiles: too far in time")

    if model_profiles.dataset.dims["id_dim"] > 0:
        ## Compute stratification metrics following observation method
        pa = coast.ProfileStratification(model_profiles)
        Zmax = 200  # metres

        pa.calc_pea(model_profiles, nemo, Zmax, CT_AS=True) #, rmax=25.0, limits=limits)
        # pa.match_to_grid(nemo,rmax=25.0,limits=limits)#,grid_name="eorca025") # this over rights nemo

        pa.dataset['delta_sss'] = pa.dataset.sss - profile.dataset.sss
        pa.dataset.delta_sss.attrs['description'] = "Sea Surface Salinity difference (model - obs)"
        pa.dataset.delta_sss.attrs['units'] = "psu"

        pa.dataset['delta_sst'] = pa.dataset.sst - profile.dataset.sst
        pa.dataset.delta_sst.attrs['description'] = "Sea Surface Temperature difference (model - obs)"
        pa.dataset.delta_sst.attrs['units'] = "deg C"

        pa.dataset['delta_pea'] = pa.dataset.pea - profile.dataset.pea
        pa.dataset.delta_pea.attrs['description'] = "Potential Energy Anomaly difference (model - obs)"
        pa.dataset.delta_pea.attrs['units'] = "J / m^3"

        fn_out = '{0}/{1}_{2}_{3}_NEMO_PEA_SST_SSS_binned.nc'.format(out_path, LME_Name, yr_str, mm_str)
        pa.dataset.to_netcdf(fn_out)

    else:
        print(f"No profiles to output")
    stop_time = time.perf_counter()
    print(f"Time elapsed: ", stop_time - start_time)
