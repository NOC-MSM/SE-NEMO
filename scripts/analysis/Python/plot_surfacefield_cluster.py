# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 16:51:43 2022

@author: Jason
"""
import numpy as np
import matplotlib.pylab as plt
import scipy.io
import pandas as pd
import re
import sys

import cartopy.crs as ccrs  # mapping plots
import cartopy.feature  # add rivers, regional boundaries etc
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER  # deg symb
from cartopy.feature import NaturalEarthFeature  # fine resolution coastline



#sys.path.insert(0, 'C:\\Users\\Jason.Goliath\\Documents\\GitHub\\COAsT\\')
sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
import coast


LME_Clusters='../Data/LME_Clusters_eORCA025.csv'
clusters = pd.read_csv(LME_Clusters)


LME_Data=np.load('../Data/LME_gridinfo_V4.npz')
a=scipy.io.loadmat('../Data/ORCA025_ROAM_GLB_LMEmaskV4.mat')
J_offset=186
LME_mask=a['LME_mask'][:,:].T
plt.figure(figsize=[11.69,8.27])
A=np.load('Position.npz')
Position=A['arr_0']
bounds=np.zeros((clusters.values.shape[0],4))
for icluster in range(clusters.values.shape[0]):
    lims=clusters.values[icluster,2:6]
    xylims=clusters.values[icluster,6:10]
    
    lims=np.array(clusters.values[icluster,2:6],dtype=int)
    Lims=np.copy(lims)
    Lims[2]=lims[2]-J_offset
    Lims[3]=lims[3]-J_offset
    
    M=LME_mask[Lims[2]:Lims[3]+1,Lims[0]:Lims[1]+1]
    #%%
    
    bathyname='../Data/eORCA025_bathy_meter.nc'
    config='example_nemo_grid_t.json'
    #bathy=coast.Gridded(fn_data=bathyname,config=config)
    #lon=bathy.dataset.longitude
    #lat=bathy.dataset.latitude
    
    
    fn_data='/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/ORCA025-SE-NEMO_1990_2019_EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2_SST_SSS_PEA_MonClimate.nc'
    fn_domain='/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/config/domain_cfg.nc'
    #nemo_clim=coast.Gridded(fn_data=fn_data,fn_domain=fn_domain,config=config)
    
    #%%
    I=np.where(~pd.isnull(clusters.values[icluster,10:]))[0]
    LMEs=[]
    for i in I:
        LMEs.append(int(re.findall(r'\d+',clusters.values[icluster,10+i])[0]))
    MM=np.ones(M.shape)*np.nan
    for LME in LMEs:
        MM[M==LME]=LME
    
#    nemo_clim_cluster=nemo_clim.subset_as_copy(y_dim=range(lims[2],lims[3]),x_dim=range(lims[0],lims[1]))
    nemo_clim_cluster=coast.Gridded(fn_data=fn_data,fn_domain=fn_domain,config=config,lims=lims)

    PEA=np.ma.masked_where(np.repeat(nemo_clim_cluster.dataset.bottom_level.values[np.newaxis,:,:]==0,12,axis=0),
                           nemo_clim_cluster.dataset.PEA_monthy_clim.values)
    
        
    #%%
#    plt.figure(icluster)#,figsize=[11.69,8.27])
    
    if True:
        x=nemo_clim_cluster.dataset.longitude.values
        y=nemo_clim_cluster.dataset.latitude.values
        central_longitude=0.0
        if np.max(x)-np.min(x)>350:
            x[x<0]=x[x<0]+360
            ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
        else:    
            ax = plt.axes(projection=ccrs.PlateCarree())
            
 
        plt.pcolormesh(x,y, PEA[7,:,:].squeeze(),transform=ccrs.PlateCarree())

        ax.set_extent(xylims,crs=ccrs.PlateCarree())
        ax.set_position(Position[icluster,:])
        #bounds[icluster ,:] =ax.get_position().bounds
        
        #plt.colorbar(orientation='vertical')
        #coastline = NaturalEarthFeature(category="physical", scale="50m", facecolor="none", name="coastline")
        #ax.add_feature(coastline, edgecolor="gray")
    
        #gl = ax.gridlines(
        #        crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.5, color="gray", alpha=0.5, linestyle="-"
        #    )
        

    plt.title('{0} {1}'.format(icluster+1,clusters.values[icluster,1]),fontsize=8)
    name=clusters.values[icluster,1].replace(' ','_')
    #plt.tight_layout()
    #plt.savefig('../Figures/'+name+'_PEA.png',bbox_inches='tight')
    #np.savez('bounds.npz',bounds)
    

    