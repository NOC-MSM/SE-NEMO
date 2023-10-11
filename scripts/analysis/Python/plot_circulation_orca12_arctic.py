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


ystart=[1991,2091]
ystop =[2000,2100]

nemo_U={}

#iexp=1
for iexp in [0,1]:# [0,2,7]:#[1,2,3,4,5,0]:
#%%    

    nemo_dir='/gws/nopw/j04/class_vol1/CLASS-MEDUSA/OUT_eORCA12/C001/monthly/'

    fn_nemo_dom='/gws/nopw/j04/class_vol1/CLASS-MEDUSA/OUT_eORCA12/C001/domain/domain_cfg.nc'
    fn_nemo_dat_t=[]
    fn_nemo_dat_u=[]
    fn_nemo_dat_v=[]
    for iy in range(ystart[iexp],ystop[iexp]+1):
        for im in range(1,12+1):    
            fn_nemo_dat_t.append(f'{nemo_dir}/{iy}/eORCA12_MED_UKESM_y{iy}m{im:02}_grid_T.nc')
            fn_nemo_dat_u.append(f'{nemo_dir}/{iy}/eORCA12_MED_UKESM_y{iy}m{im:02}_grid_U.nc')
            fn_nemo_dat_v.append(f'{nemo_dir}/{iy}/eORCA12_MED_UKESM_y{iy}m{im:02}_grid_V.nc')
        
    
    #nemo_dom= coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_t_grid,multiple=False)
    nemo_t=coast.CurrentsOnT(fn_domain=fn_nemo_dom,config='example_nemo_grid_t.json',multiple =True,no_depths=True)
    nemo_u = coast.Gridded(fn_data=fn_nemo_dat_u, config=fn_config_u_grid,multiple=True)
    nemo_v = coast.Gridded(fn_data=fn_nemo_dat_v, config=fn_config_v_grid,multiple=True)
    nemo_t.subset(z_dim=[0])
    nemo_u.subset(z_dim=[0])
    nemo_v.subset(z_dim=[0])

    #nemo_t.currents_on_t(nemo_u,nemo_v)
    #NNA
#   x_min=-79;x_max=12;y_min=26;y_max=69
#Atlantic
    x_min=-45;x_max=31;y_min=55;y_max=90
    
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
    Name='{2}-{0}-{1}'.format(ystart[iexp],ystop[iexp],'ORCA0083')
    
    #plotting
    SP,US,VS=surf.mean_surface_circulation(nemo_u,nemo_v,nemo_t,mask)
    surf.plot_surface_circulation(SP,US,VS,nemo_t,mask,Name,Vmax=0.3,
                                           Np=12
                                           ,headwidth=5,scale=80,minshaft=2      
                                           )
    plt.savefig('../Figures/Circulation/Surface_Currents_Arctic_' + Name.replace(' ','_')+'.png',dpi=300)
    #fn_out=("/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/Circulation/Surface_Currents_NEA_{0}.nc".format(Name)).replace(' ','_')
    #nemo_t_out=coast.Gridded(nemo_t,config=fn_config_t_grid)
    #surf.save_currents(SP,US,VS,fn_out,nemo_t_out)
    nemo_U[iexp]=nemo_u
    
    
    
