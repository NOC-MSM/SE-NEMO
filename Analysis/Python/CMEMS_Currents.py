#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 12:28:05 2022

@author: jholt
"""
import socket
isliv = 'livljobs' in socket.gethostname()
if isliv:
    coast_dir='/login/jholt/work/Git/COAsT/'
else:
    coast_dir='/home/users/jholt/work/Git/COAsT/'
import sys
sys.path.insert(0,coast_dir)
import coast
import numpy as np
from getpass import getpass
import circulation
USERNAME='jholt'
PASSWORD=getpass( 'Password: ' )

database = coast.Copernicus(USERNAME, PASSWORD, "my")
globcurrent=database.get_product("dataset-uv-rep-monthly")

nemo_t = coast.Gridded(fn_data=globcurrent, config="/home/users/jholt/work/Git/COAsT/config/example_cmems_grid_uv.json")
nt=nemo_t.dataset.t_dim.size
T=np.array([])
for im in [8,9,10]:
    T=np.append(T,np.arange(im,nt,12))
T=np.sort(T).astype(int)    

nemo_t.subset(x_dim=range(630,770),y_dim=range(520,630),z_dim=0,t_dim=T)
mask=nemo_t.dataset.u_velocity[0,:,:].values != np.nan
circulation.plot_surface_circulation(nemo_t, nemo_t,nemo_t, mask,'CMEMS', co_located=True,Vmax=.32)