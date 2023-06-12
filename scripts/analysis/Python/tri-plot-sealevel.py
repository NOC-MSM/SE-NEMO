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

import matplotlib.pylab as plt
import coast
import numpy as np
import pickle
import pandas as pd
import cartopy.crs as ccrs  # mapping plots
import surfacefields as sf 
cmap1=sf.lightcolormap(32,2)
cmap1.set_bad([0.75,0.75,0.75])



names,dpaths,DOMS,_  = coast.experiments(experiments='experiments-CLASS-triad.json')

year_start=1990
year_stop =1990
ssh_mean={}
x={}
y={}
for iexp in [0,1,2]:
    
    grid='T'
    
    january = 1
    december = 13  # range is non-inclusive so we need 12 + 1
    x_min=-85
    x_max=13
    y_min=26
    y_max=70

    fn_config_t_grid='../Config/senemo_grid_t.json'
    directory=dpaths[iexp]
    run_name=names[iexp]
    nemo_dom=coast.Gridded(fn_domain= DOMS[iexp],config=fn_config_t_grid,no_depths=True)
    j,i,_=nemo_dom.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
    imin=min(i)
    imax=max(i)
    jmin=min(j)
    jmax=max(j)
    lims=[imin,imax,jmin,jmax] 
    fnames = []
    for year in range(year_start, year_stop + 1):
            for month in range(january, december):
                new_name = f"{directory}/{year}/e{run_name}_MED_UKESM_y{year}m{month:02}_grid_{grid}.nc"
                fnames.append(new_name)
    fnames=np.array(fnames)
    
    nemo = coast.Gridded(fn_data= fnames, fn_domain= DOMS[iexp],config=fn_config_t_grid,multiple=True,lims=lims)
    
    #nemo.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    ssh_mean[iexp]=np.ma.masked_where(nemo.dataset.bottom_level==0,
                                      np.mean(nemo.dataset.ssh.values[:,:,:],axis=0))
    x[iexp]=nemo.dataset.longitude.values
    y[iexp]=nemo.dataset.latitude.values
    
        #%%
fig, axs = plt.subplots(nrows=3,ncols=1,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(8.5,11))

for iexp in [0,1,2]:
    X=x[iexp]
    Y=y[iexp]
    VAR=ssh_mean[iexp]
    vmin=-1
    vmax=1
    xylims=[x_min,x_max,y_min,y_max]
    
    im=axs[iexp].pcolormesh(X,Y, VAR,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap1)
    axs[iexp].set_extent(xylims,crs=ccrs.PlateCarree())
    #plt.colorbar(im,orientation='horizontal')
    axs[iexp].set_title(names[iexp])

cbar_ax = fig.add_axes([0.85, 0.11, 0.05, 0.25])

# Draw the colorbar
cbar=fig.colorbar(im, cax=cbar_ax,orientation='vertical')

