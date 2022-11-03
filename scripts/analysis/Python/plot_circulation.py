#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 16:33:37 2022

@author: jholt
"""

coast_dir='/home/users/jholt/work/Git/COAsT/'
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

nemo_dir='/gws/nopw/j04/class_vol1/CLASS-MEDUSA/OUT_eORCA12/C001/monthly/'
fn_nemo_dat_u= coast.nemo_filenames(nemo_dir,'CLASS-ORCA0083',ystart,ystop,grid='U') 
fn_nemo_dat_v= coast.nemo_filenames(nemo_dir,'CLASS-ORCA0083',ystart,ystop,grid='V') 
#%%
fn_nemo_dom='/gws/nopw/j04/class_vol1/CLASS-MEDUSA/OUT_eORCA12/C001/domain/domain_cfg.nc'
fn_nemo_dat_t=nemo_dir+'2012/eORCA12_MED_UKESM_y2012m12_grid_T.nc'


nemo_dom= coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_t_grid,multiple=False)
nemo_t = coast.Gridded(fn_data=fn_nemo_dat_t, config=fn_config_t_grid,multiple=True)
nemo_u = coast.Gridded(fn_data=fn_nemo_dat_u, config=fn_config_u_grid,multiple=True)
nemo_v = coast.Gridded(fn_data=fn_nemo_dat_v, config=fn_config_v_grid,multiple=True)

#%%
x_min=-28;x_max=12;y_min=40;y_max=69

j,i,_=nemo_t.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
imin=min(i)
imax=max(i)
jmin=min(j)
jmax=max(j) 
nemo_t.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
nemo_u.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
nemo_v.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
mask=(nemo_t.dataset.salinity>0).squeeze()
Name='ORCA12-{0}-{1}'.format(ystart,ystop)
#%%
#plotting
SP,US,VS=circ.mean_surface_circulation(nemo_u,nemo_v,nemo_t,mask)
surf.plot_surface_circulation(SP,US,VS,nemo_t,mask,Name,
                                       Np=6
                                       ,headwidth=5,scale=80,minshaft=2      
                                       )
plt.savefig('../Figures/Circulation/Surface_Currents_NEA_' + Name.replace(' ','_')+'.png',dpi=300)
fn_out=("/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/Circulation/Surface_Currents_NEA_{0}.nc".format(Name)).replace(' ','_')
nemo_t_out=coast.Gridded(nemo_t,config=fn_config_t_grid)
circ.save_currents(SP,US,VS,fn_out,nemo_t_out)
