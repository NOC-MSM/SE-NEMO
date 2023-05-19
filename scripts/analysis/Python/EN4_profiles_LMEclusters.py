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

import socket
isliv = 'livljobs' in socket.gethostname()
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
import coast
LME_Clusters='../Data/LME_Clusters_eORCA025.csv'
clusters = pd.read_csv(LME_Clusters)
nclusters=clusters.values.shape[0]
fn_profile_config='../Config/example_en4_profiles.json'
en4_path='/projectsa/NEMO/OBS/EN4.2.1/'
en4_name='EN.4.2.1.f.profiles.l09.'
yr=np.ones(1)*np.nan
#Years = np.arange(1970,2018+1)
Years = np.arange(1980,2018+1)
Years=[1980]
Months =np.arange(1,12+1)

fn_profile=[]
nemo=coast.Gridded(fn_domain='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc',config='example_nemo_grid_t.json')
for iy in Years:
  for im in Months:   

   yy=str(iy)
   mm=str(im)
   if im<10:
    mm='0'+mm   
   fn_profile.append(en4_path+en4_name+yy+mm+'.nc')
fn_profilenames=np.array(fn_profile)   
for icluster in [11]:#range(nclusters):
    print(clusters.values[icluster,1])
    lonlims=clusters.values[icluster,6:8]
    latlims=clusters.values[icluster,8:10]
    
    profile = coast.Profile(config=fn_profile_config)
    profile.read_en4(fn_profile,multiple=True)
    profile_sub=profile.subset_indices_lonlat_box(lonlims,latlims)
    profile_sub.dataset = profile.dataset.isel(id_dim=np.arange(0, profile_sub.dataset.dims["id_dim"], 100)).load()
    pa = coast.ProfileStratification(profile_sub)
    Zmax = 200  # metres
    pa.clean_data(profile_sub,nemo,Zmax)
    pa.calc_pea(profile_sub,nemo, Zmax)
    fname='../Data/{0}_EN4_PEA_SST_SSS_test.nc'.format(clusters.values[icluster,1].replace(' ','_'))
    pa.to_netcdf(fname)
    