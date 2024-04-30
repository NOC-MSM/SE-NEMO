#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Companion script to find_nemo_profile_by_LME.ipynb (which should be called find_EN4_profiles_on_nemo_grid_by_LME.py)

For a single LME extracts synthetic EN4 profiles from NEMO, iterating over monthly files and concatenating over id_dim
Saves SST, SSS, PEA and difference from EN4

Master EN4 profiles already extracted by
EN4_profiles_LME.py

Created on Mar 20 2024

target machine: JASMIN / lotus
@author: jelt

ssh sci6 -L 6000:localhost:22
pycharm jelt@loclahost:6000
env: senemo-profile-lme (over ssh)


#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
rm LOGS/OUT* LOGS/*.err LOGS/*.out
for (( iLME=34; iLME<36; iLME++ ))
do
 echo "sbatch $iLME lotus_NEMO_profiles_by_LME.sh $iLME"
 sbatch -J $iLME lotus_NEMO_profiles_by_LME.sh $iLME
done


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


def extract_profiles_from_gridded(profile, gridded):
    """
    input: profile object containing indices where located in the to gridded object
            gridded object containing target data
    output: extracted_gridded_profiles, profile_matched - new synthetic profiles and input profiles, trimmed over time to match
    """
    extracted = coast.Profile(config=fn_profile_config)
    varlist = ["conservative_temperature", "absolute_salinity", "depth", "e3_0"]
    store2d = np.ones((len(varlist), profile.dataset.dims["id_dim"], gridded.dataset.dims["z_dim"] ))
    store1d = np.ones((3, profile.dataset.dims["id_dim"]))
    time_max = gridded.dataset.time.max().values
    time_min = gridded.dataset.time.min().values
    if time_max == time_min: # stretch time envelope to catch monthly EN4 data mapped to around the middle of the month
        time_max = time_max + np.timedelta64(5, "D")
        time_min = time_min - np.timedelta64(5, "D")
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
            #print('out of time range - skip')

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



# read in specified iLME
args = sys.argv
try:
    iLME = int(args[1])
except:
    iLME = 34  # Gulf of Thailand
    print('No input LME specified. Use Gulf of Thailand')

sys.path.insert(0,'/home/users/jelt/GitHub/COAsT/')

import coast



ystart = 1978
ystop = 2019

# Load the rebinned EN4 data by LME. These are averaged onto the NEMO space-time grid
LME_en4_in_dir = "/gws/nopw/j04/class_vol2/senemo/jelt/PROCESSED/EN4.2.1/1978-2019/"
in_path = LME_en4_in_dir
out_path = LME_en4_in_dir

#%% load processed EN4 data
# load LME data for names
LME_Data=np.load('/home/users/jelt/GitHub/SE-NEMO/scripts/analysis/Data/LME_gridinfo_V4.npz')

nLME = LME_Data["DOMNAM"].shape[0]
#for iLME in range(34,36): #66):
LME_Name=LME_Data["DOMNAM"][iLME]
print(f'{LME_Name}: {iLME} of {nLME}')

start_time = time.perf_counter()
first_time = True

fn_in = '{0}/{1}_{2}_{3}_EN4_PEA_SST_SSS_binned.nc'.format(in_path, LME_Name, ystart, ystop)
fn_profile_config='/home/users/jelt/GitHub/COAsT/config/example_en4_profiles.json'

profile = coast.Profile(config=fn_profile_config)
profile.dataset = xr.open_dataset(fn_in, chunks={'id_dim': 10000})

## Do we need z_dim ?
#profile.dataset = profile.dataset.expand_dims(dim={"z_dim": 1})
#profile.dataset['depth'] = xr.DataArray(np.zeros((profile.dataset.dims["id_dim"], 1)), dims = ("id_dim", "z_dim"))


# Set up stationary NEMO paths
#fn_nemo_data = "/gws/nopw/j04/class_vol2/senemo/FINAL_PHYSICS_RUNS/GS1p0_notide/output/SENEMO_1m*grid_T*.nc"
fn_nemo_grid = "/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/config/domain_cfg.nc"
fn_nemo_config_t = "/gws/nopw/j04/class_vol2/senemo/jelt/PROCESSED/tmp_nemo_grid_t.json"

## Loop over year and load NEMO data if there are profiles to extract
#for year in range(1985, 1986+1):
for year in range(ystart, ystop+1):
    for month in range(1, 12+1):
        yr_mon_str = str(year) + '-' + str(month).zfill(2)
        #if (year in profile.dataset.time.dt.year) and (month in profile.dataset.time.dt.month):
        if np.datetime64(yr_mon_str) in profile.dataset.time.astype('datetime64[M]'):
            fn_nemo_data = "/gws/nopw/j04/class_vol2/senemo/FINAL_PHYSICS_RUNS/GS1p0_notide/output/SENEMO_1m_*grid_T*"+str(year)+str(month).zfill(2)+".nc"

            print(f'load NEMO month-year: {yr_mon_str}')
            print(f"event matches: {(np.datetime64(yr_mon_str) == profile.dataset.time.astype('datetime64[M]').values).sum()}")
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


            # For the global model the wrapped seam can cause problems with the nearest neighbour calculation if the NEMO seam is within the LME
            # if ((nemo.dataset.longitude.isel(x_dim=0).mean() < profile.dataset.longitude.values.max())
            #    & (nemo.dataset.longitude.isel(x_dim=0).mean() > profile.dataset.longitude.values.min())):
            #    print(f"Need to roll the dataset {LME_Name}")
            #    nemo.dataset = nemo.dataset.roll(x_dim=nemo.dataset.dims["x_dim"]//2, roll_coords=True)

            # Use obs operator for horizontal remapping of Gridded onto Profile.
            if(0):
                model_profiles = profile.obs_operator(nemo)

                ### Discard profiles where the interpolation distance is too large
                too_far = 13.75  # 110/4/2 distance km
                keep_indices = model_profiles.dataset.interp_dist <= too_far
                model_profiles = model_profiles.isel(id_dim=keep_indices)
                profile = profile.isel(id_dim=keep_indices)
                if np.sum(~keep_indices.values) > 0:
                    print(f"Dropped {np.sum(~keep_indices.values)} of {len(keep_indices)} profiles: too far in space")

                # Throw away profile where the interpolation time is larger than 48h
                keep_indices = np.abs(model_profiles.dataset.interp_lag) <= np.timedelta64(48,'h')  ## SHOULD MOVE PARAMTER TO CONFIG
                model_profiles = model_profiles.isel(id_dim=keep_indices)
                profile = profile.isel(id_dim=keep_indices)
                if np.sum(~keep_indices.values) > 0:
                    print(f"Dropped {np.sum(~keep_indices.values)} profiles: too far in time")

            else:
                model_profiles, subset_profile = extract_profiles_from_gridded(profile, nemo)


            if model_profiles.dataset.dims["id_dim"] > 0:
                ## Compute stratification metrics following observation method
                pa = coast.ProfileStratification(model_profiles)
                Zmax = 200  # metres

                pa.calc_pea(model_profiles, nemo, Zmax, CT_AS=True) #, rmax=25.0, limits=limits)
                # pa.match_to_grid(nemo,rmax=25.0,limits=limits)#,grid_name="eorca025") # this over rights nemo

                pa.dataset['delta_sss'] = pa.dataset.sss - subset_profile.dataset.sss
                pa.dataset.delta_sss.attrs['description'] = "Sea Surface Salinity difference (model - obs)"
                pa.dataset.delta_sss.attrs['units'] = "psu"

                pa.dataset['delta_sst'] = pa.dataset.sst - subset_profile.dataset.sst
                pa.dataset.delta_sst.attrs['description'] = "Sea Surface Temperature difference (model - obs)"
                pa.dataset.delta_sst.attrs['units'] = "deg C"

                pa.dataset['delta_pea'] = pa.dataset.pea - subset_profile.dataset.pea
                pa.dataset.delta_pea.attrs['description'] = "Potential Energy Anomaly difference (model - obs)"
                pa.dataset.delta_pea.attrs['units'] = "J / m^3"

                ## Create new dataset or append.
                print(f"Saving {model_profiles.dataset.dims['id_dim']} profiles")
                if first_time:
                    first_time = False
                    datasets = []
                datasets.append(pa.dataset)

            else:
                print(f"No profiles to output")

        else:
            print(f'Skipping month-year:{yr_mon_str}')


# concatenate at the end - to make it more efficient
pa_datasets = xr.concat(datasets, dim="id_dim")

# Save to file
fn_out = '{0}/{1}_{2}_{3}_NEMO_PEA_SST_SSS_binned.nc'.format(out_path, LME_Name, ystart, ystop)
print(f"Saving combined model and obs data into a single file: {fn_out}")
pa_datasets.to_netcdf(fn_out)

stop_time = time.perf_counter()
print(f"Time elapsed: ", stop_time - start_time)







