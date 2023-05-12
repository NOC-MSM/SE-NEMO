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
    coast_dir='/home/users/jholt/Git/COAsT/'
import sys
sys.path.insert(0,coast_dir)
import coast
import numpy as np
from getpass import getpass
import surfacefields
import scipy.io
import matplotlib.pylab as plt
USERNAME='jholt'
PASSWORD=getpass( 'Password: ' )

database = coast.Copernicus(USERNAME, PASSWORD, "my")
#Select global model 1/12
globcurrent=database.get_product("cmems_mod_glo_phy_my_0.083_P1M-m")
#Select Altimeter observations
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
#set time period
T=np.arange(nt)
T=np.arange(120)
T=np.arange(12)
A=np.load('../Data/LME_gridinfo_equ025.npz')
a=scipy.io.loadmat('../Data/equalgrid_025_LMEmask.mat')
nlme=66

#LME_mask=a['LME_mask'][:,:].T
lmelist=np.array([34])-1
#lmelist=np.arange(nlme)
#for ilme in lmelist:    
#        LMENAM=A['DOMNAM'][ilme]
#        x_min=A['x_min'][ilme]
#        x_max=A['x_max'][ilme]
#        y_min=A['y_min'][ilme]
#        y_max=A['y_max'][ilme]        
        
x_min=-15
x_max=13
y_min=45
y_max=65

x_min=-79
x_max=12
y_min=26
y_max=69
nemo_t.make_lonLat_2d()
j,i,_=nemo_t.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])

imin=min(i)
imax=max(i)
jmin=min(j)
jmax=max(j)        

LMENAM='NWS'
name="CMEMS Glob Current"
name="CMEMS ORCA12 Current"

SEASON='JAS'
YEARS="1993"#-2002"
nemo_t1=nemo_t.subset_as_copy(x_dim=range(imin,imax),y_dim=range(jmin,jmax),z_dim=0,t_dim=T)
mask=nemo_t1.dataset.u_velocity[0,:,:].values != np.nan

REGION=LMENAM
#REGION='NWS'  
Name=name+' '+SEASON+' '+YEARS+' '+REGION
#%%
SP,US,VS=surfacefields.mean_surface_circulation(nemo_t1,nemo_t1,nemo_t1,mask,co_located=True)         
surfacefields.plot_surface_circulation(SP, US, VS,nemo_t1, mask,Name,Vmax=.16,Np=12)

plt.savefig('../Figures/Circulation/Surface_Currents_' + Name.replace(' ','_')+'.png',dpi=300)

#fn_out=("/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/Circulation/Surface_Currents_{0}.nc".format(Name)).replace(' ','_')
#nemo_t_out=nemo_t1.copy()
#for var in list(nemo_t_out.dataset.keys()):
#    print(var)
#    nemo_t_out.dataset=nemo_t_out.dataset.drop(var)  
 
#circulation.save_currents(SP,US,VS,fn_out,nemo_t_out)        
