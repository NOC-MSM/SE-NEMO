#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 15:56:54 2023

@author: jholt
"""
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import sys
sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
import coast
import cartopy.crs as ccrs  # mapping plots

config='example_nemo_grid_t.json'
fn_domain='/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/config/domain_cfg.nc'
LME_Clusters='../Data/LME_Clusters_eORCA025.csv'
clusters = pd.read_csv(LME_Clusters)
nclusters=clusters.values.shape[0]

bounds=np.zeros((nclusters,4))
Dlat=np.zeros(nclusters)
for icluster in range(nclusters):
    plt.figure(0)
    
    lims=clusters.values[icluster,2:6]
    xylims=clusters.values[icluster,6:10]
    
    lims=np.array(clusters.values[icluster,2:6],dtype=int)
    Lims=np.copy(lims)
    #Lims[2]=lims[2]-J_offset
    #Lims[3]=lims[3]-J_offset
    

    nemo_domain=coast.Gridded(fn_domain=fn_domain,config=config,lims=lims)
    D=nemo_domain.dataset.bathymetry.values.squeeze()
    x=nemo_domain.dataset.longitude.values
    y=nemo_domain.dataset.latitude.values
    central_longitude=0.0
    if np.max(x)-np.min(x)>350:
        x[x<0]=x[x<0]+360
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    else:    
        ax = plt.axes(projection=ccrs.PlateCarree())
        
     
    plt.pcolormesh(x,y, D.squeeze(),transform=ccrs.PlateCarree())
    
    ax.set_extent(xylims,crs=ccrs.PlateCarree())
    gl = ax.gridlines(
                crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.5, color="gray", alpha=0.5, linestyle="-"
            )
    bounds[icluster ,:] =ax.get_position().bounds
    Dlat[icluster]=xylims[3]-xylims[2]
    plt.close()
#%%
SC=1.
mDlat=np.mean(Dlat)
Row={}
Row[1]=np.array([6,5,4,2,3]).astype(int)-1
Row[2]=np.array([14,23,15,22,1]).astype(int)-1
Row[3]=np.array([11,16,21,17,7]).astype(int)-1
Row[4]=np.array([12,20,18,19,8,9,10]).astype(int)-1
Row[5]=np.array([13]).astype(int)-1

plt.figure(figsize=[11.69,8.27])
gap=0.05
Position=np.zeros((23,4))

for ir in range(1,6):
    #A=np.load('bounds.npz')
    #bounds=A['arr_0']

    bnds=bounds[Row[ir],:]
    nP=np.size(Row[ir])
    sc=Dlat[Row[ir]]/np.mean(Dlat)
    sc=1
    if ir != 5:
        WW=np.sum(sc*bnds[:,2])+(nP-1)*gap+2*gap
    gap_s=gap/WW

    x0=gap_s
    if ir ==1 :
        y0 = 1-x0*1.5

    Hmax=0.
    for IP,ip in enumerate(Row[ir]):
#        W=SC*sc[IP]*bounds[ip,2]/WW
        W=bounds[ip,2]/WW
        H=W*bounds[ip,3]/bounds[ip,2]
        y1=y0-H
        Position[ip,0] = x0
        Position[ip,1] = y1
        Position[ip,2] = W
        Position[ip,3] = H    
        x0=x0+gap_s+W
        ax=plt.axes()
        ax.set_position(Position[ip,:])
        Hmax=np.max([H,Hmax])
    print(y0)    
    y0=y0 - Hmax -gap_s*1.5
np.savez('Position.npz',Position)
