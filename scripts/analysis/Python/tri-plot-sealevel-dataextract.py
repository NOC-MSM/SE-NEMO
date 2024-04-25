#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 16:19:06 2023

@author: jholt
"""
import sys
import socket
isliv = 'livljobs' in socket.gethostname()
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
sys.path.insert(0,'./')
import matplotlib.pylab as plt
import coast
import numpy as np
import pickle
import pandas as pd
import cartopy.crs as ccrs  # mapping plots
import surfacefields as sf
import scipy
import xarray as xr
cmap1=sf.lightcolormap(32,2)
cmap1.set_bad([0.75,0.75,0.75])

x_min = -85
x_max = 13
y_min = 26
y_max = 70

x_min = -19
x_max = 13
y_min = 40
y_max = 65

names,dpaths,DOMS,_  = coast.experiments(experiments='experiments-CLASS-triad2.json')

#%%
year_start=1990
year_stop =2009
ssh_mean={}
ssh_std={}

x={}
y={}
for iexp in [2]:#5,4,3]:#[0,1,2,5,4,3]:#[0,1,2]:
    
    grid='T'
    
    january = 1
    december = 13  # range is non-inclusive so we need 12 + 1


    fn_config_t_grid='../Config/senemo_grid_t.json'
    directory=dpaths[iexp]
    run_name=names[iexp]
    nemo_dom=coast.Gridded(fn_domain= DOMS[iexp],config=fn_config_t_grid,multiple=False,no_depths=True)
    j,i,_=nemo_dom.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
    imin=min(i)
    imax=max(i)
    jmin=min(j)
    jmax=max(j)
    lims=[imin,imax,jmin,jmax] 
    fnames = []

    if iexp <3 :
        for year in range(year_start, year_stop + 1):
                for month in range(january, december):
                    if iexp==0:
                        new_name = f"{directory}/{year}/{run_name}_N06_{year}{month:02}m01T.nc"
                    else:
                        new_name = f"{directory}/{year}/{run_name}-N06_{year}m{month:02}T.nc"
                    fnames.append(new_name)
        fnames=np.array(fnames)
    else:
        fnames= coast.nemo_filename_maker(dpaths[iexp],year_start,year_stop,grid='T') 
    nemo = coast.Gridded(fn_data= fnames, fn_domain= DOMS[iexp],config=fn_config_t_grid,multiple=True,lims=lims,no_depths=True)
    
    #nemo.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    #ssh_mean[iexp]=np.ma.masked_where(nemo.dataset.bottom_level==0,
    #                                  np.mean(nemo.dataset.ssh.values[:,:,:],axis=0))
    #ssh_std[iexp]=np.ma.masked_where(nemo.dataset.bottom_level==0,
    #                                  np.std(nemo.dataset.ssh.values[:,:,:],axis=0))

    ll=list(nemo.dataset.keys())
    ll.remove('ssh')
    ll.remove('time_centered')
    ll.remove('bottom_level')
    dataset_new=nemo.dataset.drop_vars(ll)
    ssh=dataset_new.ssh.values
    SSH=xr.DataArray(ssh,dims=dataset_new.dims,coords=dataset_new.coords)
    dataset_new2=dataset_new.drop_vars('ssh')
    dataset_new2['ssh']=SSH
    outname=f"/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/{names[iexp]}_ssh_1990-2009.nc"
    dataset_new2.to_netcdf(outname)
