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
import scipy.io

USERNAME='jholt'
PASSWORD=getpass( 'Password: ' )

database = coast.Copernicus(USERNAME, PASSWORD, "my")
globcurrent=database.get_product("cmems_mod_glo_phy_my_0.083_P1M-m")
#globcurrent=database.get_product("dataset-uv-rep-monthly")


#%%
nemo_t = coast.Gridded(fn_data=globcurrent, 
                       config="example_cmems_grid_uv.json",
                       Make_LonLat_2D=True)
nt=nemo_t.dataset.t_dim.size
#nt=12
T=np.array([])

#select months
for im in [6,7,8]:
    T=np.append(T,np.arange(im,nt,12))
T=np.sort(T).astype(int)    

A=np.load('../Data/LME_gridinfo_equ025.npz')
a=scipy.io.loadmat('../Data/equalgrid_025_LMEmask.mat')
nlme=66

LME_mask=a['LME_mask'][:,:].T
lmelist=np.array([34])-1
for ilme in lmelist:    

        x_min=A['x_min'][ilme]
        x_max=A['x_max'][ilme]
        y_min=A['y_min'][ilme]
        y_max=A['y_max'][ilme]        
                
#        x_min=-15
#        x_max=13
#        y_min=45
#        y_max=65
        
        
        j,i,_=nemo_t.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
        imin=min(i)
        imax=max(i)
        jmin=min(j)
        jmax=max(j)        
#        jmin=375
#        jmax=455
#        imin=1030
#        imax=1130

#        imin=1950
#        imax=2300
#        jmin=1450
#        jmax=1750
#        LMENAM='_ORCA12'
        nemo_t1=nemo_t.subset_as_copy(x_dim=range(imin,imax),y_dim=range(jmin,jmax),z_dim=0,t_dim=T)
        mask=nemo_t1.dataset.u_velocity[0,:,:].values != np.nan
        circulation.plot_surface_circulation(nemo_t1, nemo_t1,nemo_t1, mask,'CMEMS ORCA12 '+LMENAM, co_located=True,Vmax=.32)
