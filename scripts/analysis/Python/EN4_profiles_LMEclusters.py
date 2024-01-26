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

import coast

coast_path='/login/jholt/work/Git/COAsT/'
en4_path='/projectsa/NEMO/OBS/EN4.2.1/'
en4_name='EN.4.2.1.f.profiles.l09.'

LME_Clusters='../Data/LME_Clusters_eORCA025.csv'
clusters = pd.read_csv(LME_Clusters)
nclusters=clusters.values.shape[0]
fn_profile_config='../Config/example_en4_profiles.json'

nemo=coast.Gridded(fn_domain='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc',config='example_nemo_grid_t.json')

ystart = 2019
ystop = 2019

for icluster in [11]:#range(nclusters): # needs to work for cluster 17
    start_time = time.perf_counter()
    print(clusters.values[icluster,1])
    lonlims=clusters.values[icluster,6:8]
    latlims=clusters.values[icluster,8:10]
    x_min,x_max=lonlims
    y_min,y_max=latlims
    region_bounds = [x_min,x_max,y_min,y_max]
    j,i,_=nemo.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
    limits=[min(j),max(j),min(i),max(i)]
    profile = coast.Profile(config=fn_profile_config)
    fn_profile = profile.make_filenames(en4_path, 'EN4', ystart, ystop)
    profile= profile.extract_en4_profiles(fn_profile, region_bounds,chunks={"id_dim" : 100})

    pa = coast.ProfileStratification(profile)
    Zmax = 200  # metres

    pa.calc_pea(profile,nemo, Zmax)
    pa.match_to_grid(nemo,limits=limits,grid_name="eorca025") # this over rights nemo
    fname='../Data/{0}_EN4_PEA_SST_SSS_test.nc'.format(clusters.values[icluster,1].replace(' ','_'))
    pa.dataset.to_netcdf(fname)
    stop_time = time.perf_counter()
    print(stop_time - start_time)
