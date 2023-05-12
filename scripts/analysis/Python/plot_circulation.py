#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 16:33:37 2022

@author: jholt
"""

coast_dir='/home/users/jholt/Git/COAsT/'
import sys
sys.path.insert(0,coast_dir)
import matplotlib.pylab as plt
import coast
import numpy as np
import xarray as xr
import scipy.io
import surfacefields as surf
fn_config_u_grid='./example_nemo_grid_u.json'
fn_config_v_grid='./example_nemo_grid_v.json'
fn_config_t_grid='./example_nemo_grid_t.json'

ystart=2090
ystop=2099
ystart=1991
ystop=2000

names,dpaths,DOMS,_  = coast.experiments(experiments='experiments.json')
#nemo_U={}

#iexp=1
for iexp in [3]:
#%%    
    print(names[iexp])
    #nemo_dir='/gws/nopw/j04/class_vol1/CLASS-MEDUSA/OUT_eORCA12/C001/monthly/'
    fn_nemo_dat_u= coast.nemo_filename_maker(dpaths[iexp],ystart,ystop,grid='U') 
    fn_nemo_dat_v= coast.nemo_filename_maker(dpaths[iexp],ystart,ystop,grid='V') 
    #fn_nemo_dom='/gws/nopw/j04/class_vol1/CLASS-MEDUSA/OUT_eORCA12/C001/domain/domain_cfg.nc'
    #fn_nemo_dat_t=nemo_dir+'2012/eORCA12_MED_UKESM_y2012m12_grid_T.nc'
    
    
    #nemo_dom= coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_t_grid,multiple=False)
    nemo_t=coast.CurrentsOnT(fn_domain=DOMS[iexp],config='example_nemo_grid_t.json',multiple =True,no_depths=True)
    nemo_u = coast.Gridded(fn_data=fn_nemo_dat_u, config=fn_config_u_grid,multiple=True)
    nemo_v = coast.Gridded(fn_data=fn_nemo_dat_v, config=fn_config_v_grid,multiple=True)
    nemo_t.subset(z_dim=[0])
    nemo_u.subset(z_dim=[0])
    nemo_v.subset(z_dim=[0])

    nemo_t.currents_on_t(nemo_u,nemo_v)
    
    x_min=-79;x_max=12;y_min=26;y_max=69
    
    j,i,_=nemo_t.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
    imin=min(i)
    imax=max(i)
    jmin=min(j)
    jmax=max(j) 
    nemo_t.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    nemo_u.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    nemo_v.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    #%%
    mask=(nemo_t.dataset.bottom_level>0).squeeze()
    Name='{2}-{0}-{1}'.format(ystart,ystop,names[iexp])
    
    #plotting
#    SP,US,VS=surf.mean_surface_circulation(nemo_u,nemo_v,nemo_t,mask)
#    surf.plot_surface_circulation(SP,US,VS,nemo_t,mask,Name,
#                                           Np=6
#                                           ,headwidth=5,scale=80,minshaft=2      
#                                           )
#    plt.savefig('../Figures/Circulation/Surface_Currents_NNA_' + Name.replace(' ','_')+'.png',dpi=300)
    #fn_out=("/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/Circulation/Surface_Currents_NEA_{0}.nc".format(Name)).replace(' ','_')
    #nemo_t_out=coast.Gridded(nemo_t,config=fn_config_t_grid)
    #surf.save_currents(SP,US,VS,fn_out,nemo_t_out)
    nemo_U[iexp]=nemo_u
    
    
    
    #%%
U_cmems=nemo_t1.subset_as_copy(x_dim=100).dataset.u_velocity.mean(dim='t_dim').values
lat_cmems=nemo_t1.subset_as_copy(x_dim=100).dataset.latitude
#%%
U1=nemo_U[1].subset_as_copy(x_dim=100).dataset.u_velocity.mean(dim='t_dim').values
U2=nemo_U[2].subset_as_copy(x_dim=100).dataset.u_velocity.mean(dim='t_dim').values
U0=nemo_U[0].subset_as_copy(x_dim=100).dataset.u_velocity.mean(dim='t_dim').values
U3=nemo_U[3].subset_as_copy(x_dim=100).dataset.u_velocity.mean(dim='t_dim').values

lat1=nemo_u.subset_as_copy(x_dim=100).dataset.latitude.values

#%%
j1=120
j2=95    
plt.plot(lat1[:j1],U0.squeeze()[:j1],lat1[:j1],
                          U1.squeeze()[:j1],lat1[:j1],U2.squeeze()[:j1],
                          lat1[:j1],U3.squeeze()[:j1],
                          lat_cmems[:j2],U_cmems[:j2],'.-')
plt.title('Eastward surface current 54W')
plt.xlabel('Latitude')
plt.ylabel('m/s')
plt.legend([names[0],names[1],names[2],names[3],'CMEMS Glob Current'])
plt.savefig('../Figures/Circulation/Surface_Currents_GS_CMEMS_i100.png',dpi=300)
