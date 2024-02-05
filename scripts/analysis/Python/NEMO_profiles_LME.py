#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Companion script to EN4_profiles_LME.py

Extracts profiles from LME regions that correspond to the profiles already extracted by
EN4_profiles_LME.py

Created on Tue Apr 25 13:41:16 2023

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

import socket
isliv = 'livljobs' in socket.gethostname()
if isliv:
 sys.path.insert(0,'/login/$USER/work/GitHub/COAsT/')
 #sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')

import coast


fn_profile_config='../Config/example_en4_profiles.json'
fn_domain='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/eORCA025_bathy_meter.nc'
in_path  = '/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT/EN4.2.1/1978-2019/'
out_path = '/projectsa/NEMO/jelt/SE-NEMO/ASSESSMENT/EN4.2.1/'
ystart = 1978
ystop = 2019

#%% load model data
root = "../../../../"
# And by defining some file paths
dn_files = root + "./example_files/"
fn_nemo_dat = path.join(dn_files, "coast_example_nemo_data.nc")
fn_nemo_dom = path.join(dn_files, "coast_example_nemo_domain.nc")
fn_nemo_config = path.join(root, "./config/example_nemo_grid_t.json")

# Create gridded object:
nemo = coast.Gridded(fn_nemo_dat, fn_nemo_dom, multiple=True, config=fn_nemo_config)

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
for iLME in range(2): #66):

    LME_Name=LME_Data["DOMNAM"][iLME]

    start_time = time.perf_counter()
    print(LME_Name)

    fn_in = '{0}/{1}_{2}_{3}_EN4_PEA_SST_SSS_v1.nc'.format(in_path, LME_Name, ystart, ystop)

    profile = coast.Profile(config=fn_profile_config)
    profile.dataset = xr.open_dataset(fn_in, chunks={'id_dim': 10000})

    # Use obs operator for horizontal remapping of Gridded onto Profile.
    model_profiles = profile.obs_operator(nemo)


    ### Discard profiles where the interpolation distance is too large
    too_far = 7  # distance km
    keep_indices = model_profiles.dataset.interp_dist <= too_far
    model_profiles = model_profiles.isel(id_dim=keep_indices)

    # Throw away profile where the interpolation time is larger than 12h
    keep_indices = np.abs(model_profiles.dataset.interp_lag) <= np.timedelta64(12,
                                                                               'h')  ## SHOULD MOVE PARAMTER TO CONFIG
    model_profiles = model_profiles.isel(id_dim=keep_indices)
    profile = profile.isel(id_dim=keep_indices)
    if np.sum(~keep_indices.values) > 0:
        print(f"Dropped {np.sum(~keep_indices.values)} profiles: too far in time")


    ## Compute stratification metrics following observation method
    pa = coast.ProfileStratification(model_profiles)
    Zmax = 200  # metres

    pa.calc_pea(model_profiles, nemo, Zmax) #, rmax=25.0, limits=limits)
    # pa.match_to_grid(nemo,rmax=25.0,limits=limits)#,grid_name="eorca025") # this over rights nemo
    fn_out = '{0}/{1}_{2}_{3}_NEMO_PEA_SST_SSS_v1.nc'.format(out_path, LME_Name, ystart, ystop)
    pa.dataset.to_netcdf(fn_out)

    stop_time = time.perf_counter()
    print(stop_time - start_time)
