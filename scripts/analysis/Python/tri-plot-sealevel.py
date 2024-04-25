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

x_min=90
x_max=132
y_min=-12.8
y_max=24.7



names,dpaths,DOMS,_  = coast.experiments(experiments='experiments-CLASS-triad2.json')

#%%
year_start=1990
year_stop =2009
ssh_mean={}
ssh_std={}

x={}
y={}
for iexp in [0,1,2,5,4,3]:#[0,1,2]:
    
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
    ssh_mean[iexp]=np.ma.masked_where(nemo.dataset.bottom_level==0,
                                      np.mean(nemo.dataset.ssh.values[:,:,:],axis=0))
    ssh_std[iexp]=np.ma.masked_where(nemo.dataset.bottom_level==0,
                                      np.std(nemo.dataset.ssh.values[:,:,:],axis=0))

    x[iexp]=nemo.dataset.longitude.values
    y[iexp]=nemo.dataset.latitude.values
    
    
    
#%%
#%%
outname='/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/SL_std_SENEMO_ORCA_SEASIA1.p'
with open(outname,'wb' ) as f:
   A={}
   A['ssh_std']=ssh_std
   A['ssh_mn']=ssh_mean
   A['x']=x
   A['y']=y
#   A['slmean']=slmean
   pickle.dump(A,f)
#%%
outname='/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/SL_std_SENEMO_ORCA_SEASIA1.p'
with open(outname,'rb' ) as f:

#   A['slmean']=slmean
   A=pickle.load(f)        
   ssh_std=A['ssh_std']
   ssh_mean=A['ssh_mn']
   x=A['x']
   y=A['y']
    

#%%
fig, axs = plt.subplots(nrows=2,ncols=3,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(11,8.5))
ir=0
ic=-1
for iexp in [1,0,2,5,4,3]:
    ic=ic+1
    if ic>2:
        ic=0
        ir=1
    X=x[iexp]
    Y=y[iexp]
    if iexp == 1 or iexp== 5:
        VAR=ssh_std[iexp]
        vmin=0
        vmax=0.1
    elif iexp in [4,3]:
        VAR=ssh_std[iexp]-ssh_std[5]
        vmin=-0.05
        vmax=0.05  

    else:
        X0=x[1].ravel()
        Y0=y[1].ravel()
        ssh_std0=ssh_std[1].ravel()
        X0=X0[~ssh_std0.mask]
        Y0=Y0[~ssh_std0.mask]
        ssh_std0=ssh_std0[~ssh_std0.mask]
        
#        interp = scipy.interpolate.NearestNDInterpolator(list(zip(Y0,X0)), ssh_std0)
        interp = scipy.interpolate.LinearNDInterpolator(list(zip(Y0,X0)), ssh_std0)

        ssh_std0_in=interp(y[iexp],x[iexp])
        VAR=ssh_std[iexp]-ssh_std0_in
        vmin=-0.05
        vmax=0.05   
    xylims=[x_min,x_max,y_min,y_max]
    
    im=axs[ir,ic].pcolormesh(X,Y, VAR,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap1)
    axs[ir,ic].set_extent(xylims,crs=ccrs.PlateCarree())
    
    #plt.colorbar(im,orientation='horizontal')
#    axs[iexp].set_title(names[iexp])
#    if iexp ==2 :
    if iexp == 1 or iexp ==5:
        delta=''
    else:
        delta="$\Delta$"
    axs[ir,ic].set_title(f"{delta}STD SL {names[iexp]}")
#    else:
#        axs[iexp].set_title(f"STD SL {names[iexp]} - {names[2]}")
    if iexp == 1:
        im1=im
    if iexp == 3:
        im2=im
#    if iexp == 2:
#       cbar_ax = fig.add_axes([0.85, 0.11, 0.05, 0.20])
#       cbar=fig.colorbar(im, cax=cbar_ax,orientation='vertical')
plt.subplots_adjust(hspace=-0.3,wspace=0.2)#left=0.05,right=0.83,bottom=0.025, top=0.95)
cbar_ax = fig.add_axes([0.03, 0.375, 0.03, 0.2285])
cbar=fig.colorbar(im1, cax=cbar_ax,orientation='vertical')  

cbar_ax = fig.add_axes([0.91, 0.375, 0.03, 0.2285])
cbar=fig.colorbar(im2, cax=cbar_ax,orientation='vertical')  
plt.savefig('../Figures/SSH_std_SENEMO_ORCA-triad-1990-2009_SEASIA2.png',dpi=300)
