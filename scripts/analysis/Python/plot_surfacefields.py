#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 12:15:05 2022

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

fn_config_t_grid='./example_nemo_grid_t.json'

ystart=2090
ystop=2099
ystart=1991
ystop=2000

nemo_dir='/gws/nopw/j04/class_vol1/CLASS-MEDUSA/OUT_eORCA12/C001/monthly/'
fn_nemo_dat_t= coast.nemo_filenames(nemo_dir,'CLASS-ORCA0083',ystart,ystop,grid='T') 

#%%
fn_nemo_dom='/gws/nopw/j04/class_vol1/CLASS-MEDUSA/OUT_eORCA12/C001/domain/domain_cfg.nc'

#nemo_dom= coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_t_grid,multiple=False)
nemo_t = coast.Gridded(fn_data=fn_nemo_dat_t, config=fn_config_t_grid,multiple=True)
#%%
x_min=-28;x_max=12;y_min=40;y_max=69

j,i,_=nemo_t.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
imin=min(i)
imax=max(i)
jmin=min(j)
jmax=max(j) 
nemo_t.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))

mask=(nemo_t.dataset.salinity>0).squeeze()
Name='ORCA12-{0}-{1}'.format(ystart,ystop)

TMP,SAL=surf.mean_surface_TS(nemo_t,mask)
mask=(SAL>0).squeeze()
SAL=np.ma.masked_where(mask==0,SAL)
SALa=SAL-np.nanmean(SAL)
surf.plot_surface_T_field(SALa,nemo_t,mask,Name,Vmin=-2,Vmax=2)
plt.savefig('../Figures/Salinity/Surface_SalAnom_NEA_' + Name.replace(' ','_')+'.png',dpi=300)




