#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 13:41:16 2023

@author: jholt
"""
import numpy as np
import matplotlib.pylab as plt
import scipy.io
import pandas as pd
import re
import sys
import time

import socket
isliv = 'livljobs' in socket.gethostname()
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
coast_path='C:\\Users\\jholt\\Git\\COAsT'

sys.path.insert(0,coast_path)
import coast
LME_Data=np.load('../Data/LME_gridinfo_V4.npz')
en4_path='/projectsa/NEMO/OBS/EN4.2.1/'
en4_name='EN.4.2.1.f.profiles.l09.'

en4_path='C:\\Users\\jholt\\OneDrive - NOC\\Documents\\Data\\EN4\\'

fn_profile_config='../Config/example_en4_profiles.json'
fn_domain='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc'
fn_domain='C:\\Users\\jholt\\OneDrive - NOC\\Documents\\Data\\SENEMO\\eORCA025_bathy_meter.nc'

ystart = 2019
ystop = 2019
nLME = LME_Data["DOMNAM"].shape[0]
for iLME in  range(67): # needs to work for cluster 17
    LME_Name=LME_Data["DOMNAM"][iLME]

    nemo = coast.Gridded(fn_data=fn_domain, config='example_nemo_grid_t.json')
    start_time = time.perf_counter()
    print(LME_Name)

    x_min = LME_Data["x_min"][iLME]
    x_max = LME_Data["x_max"][iLME]
    y_min = LME_Data["y_min"][iLME]
    y_max = LME_Data["y_max"][iLME]

    region_bounds = [x_min,x_max,y_min,y_max]
    j,i,_=nemo.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
    if i[1] > i[0]:
        limits=[min(j),max(j),min(i),max(i)]
    else:
        limits = [0,0,0,0] # can't subset across longitude wrap
    profile = coast.Profile(config=fn_profile_config)
    fn_profile = profile.make_filenames(en4_path, 'EN4', ystart, ystop)
    profile= profile.extract_en4_profiles(fn_profile, region_bounds,chunks={"id_dim" : 100})

    pa = coast.ProfileStratification(profile)
    Zmax = 200  # metres

    pa.calc_pea(profile,nemo, Zmax,rmax=25.0,limits=limits)
    pa.match_to_grid(nemo,rmax=25.0,limits=limits)#,grid_name="eorca025") # this over rights nemo

    fname='../Data/{0}_EN4_PEA_SST_SSS_v1.nc'.format(LME_Name)
    pa.dataset.to_netcdf(fname)
    stop_time = time.perf_counter()
    print(stop_time - start_time)
