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
import surfacefields as sf

LME_Clusters='../Data/LME_Clusters_eORCA025.csv'
clusters = pd.read_csv(LME_Clusters)

def cluster_plot(x,y,var,vmin,vmax,Title,Figname,iexp=None):
    
    cmap1=sf.lightcolormap(32,2)
    cmap1.set_bad([0.75,0.75,0.75])
    
    plt.figure(figsize=[11.69,8.27])
    
    A=np.load('Position.npz')
    Position=A['arr_0']
        
    for icluster in range(clusters.values.shape[0]):  
        xylims=clusters.values[icluster,6:10]
        X=np.copy(x[icluster])
        Y=np.copy(y[icluster])        
        if iexp==None:
            VAR=var[icluster]
        else:
            VAR=var[icluster,iexp]
       # central_longitude=0.0    
        if np.max(X)-np.min(X)>350:
            X[X<0]=X[X<0]+360
            ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
        else:    
            ax = plt.axes(projection=ccrs.PlateCarree())
            
    
        im=plt.pcolormesh(X,Y, VAR,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap=cmap1)
    
        ax.set_extent(xylims,crs=ccrs.PlateCarree())
        ax.set_position(Position[icluster,:])
    
        
        #coastline = NaturalEarthFeature(category="physical", scale="50m", facecolor="none", name="coastline")
        #ax.add_feature(coastline, edgecolor="gray")
    
        #gl = ax.gridlines(
        #        crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.5, color="gray", alpha=0.5, linestyle=":"
        #    )
                
        plt.title('{0} {1}'.format(icluster+1,clusters.values[icluster,1]),fontsize=9,pad=-4)
        #name=clusters.values[icluster,1].replace(' ','_')
        
    
    position_ax=[0.497+0.12+0.05,0.04253,.25,0.02]
    cax=plt.axes(position=position_ax)
    plt.colorbar(im,orientation='horizontal',cax=cax)
    plt.title(Title,fontsize=10)    
    plt.savefig(Figname,bbox_inches='tight')
        #np.savez('bounds.npz',bounds)
#%%    
if __name__ == '__main__':
    LME_Data=np.load('../Data/LME_gridinfo_V4.npz')
    a=scipy.io.loadmat('../Data/ORCA025_ROAM_GLB_LMEmaskV4.mat')
    J_offset=186
    LME_mask=a['LME_mask'][:,:].T
    
    bounds=np.zeros((clusters.values.shape[0],4))
    vmax=800
    vmin=0.
    #EXP_NAM='EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2'
    #EXP_NAM='ZPS_NOTIDE'
    #%%
    x={}
    y={}
    PEA_max={}
    PEA_ann={}
    SAL_mean={}
    names,dpaths,DOMS,_  = coast.experiments(experiments='../Python/experiments.json')
    for iexp in [0,1,2]:
        for icluster in range(clusters.values.shape[0]):
            print(iexp,icluster)
        #%%
            lims=clusters.values[icluster,2:6]
        
            
            lims=np.array(clusters.values[icluster,2:6],dtype=int)
            Lims=np.copy(lims)
            Lims[2]=lims[2]-J_offset
            Lims[3]=lims[3]-J_offset
            
            M=LME_mask[Lims[2]:Lims[3]+1,Lims[0]:Lims[1]+1]
                
            bathyname='../Data/eORCA025_bathy_meter.nc'
            config='example_nemo_grid_t.json'
            #bathy=coast.Gridded(fn_data=bathyname,config=config)
            #lon=bathy.dataset.longitude
            #lat=bathy.dataset.latitude
            
            EXP_NAM=names[iexp]
            
            fn_data='/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/ORCA025-SE-NEMO_1990_2019_'+EXP_NAM+'_SST_SSS_PEA_MonClimate.nc'
            fn_domain='/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/config/domain_cfg.nc'
            #nemo_clim=coast.Gridded(fn_data=fn_data,fn_domain=fn_domain,config=config)
            
            
            I=np.where(~pd.isnull(clusters.values[icluster,10:]))[0]
            LMEs=[]
            for i in I:
                LMEs.append(int(re.findall(r'\d+',clusters.values[icluster,10+i])[0]))
            MM=np.ones(M.shape)*np.nan
            for LME in LMEs:
                MM[M==LME]=LME
            
            nemo_clim_cluster=coast.Gridded(fn_data=fn_data,fn_domain=fn_domain,config=config,lims=lims)
        
            PEA=np.ma.masked_where(np.repeat(nemo_clim_cluster.dataset.bottom_level.values[np.newaxis,:,:]==0,12,axis=0),
                                   nemo_clim_cluster.dataset.PEA_monthy_clim.values)
            
            SAL=np.ma.masked_where(np.repeat(nemo_clim_cluster.dataset.bottom_level.values[np.newaxis,:,:]==0,12,axis=0),
                                   nemo_clim_cluster.dataset.SSS_monthy_clim.values)
            
                
        
            
            x[icluster]=nemo_clim_cluster.dataset.longitude.values
            y[icluster]=nemo_clim_cluster.dataset.latitude.values
            PEA_max[icluster,iexp]=np.max(PEA,axis=0).squeeze()
            PEA_ann[icluster,iexp]=np.max(PEA,axis=0).squeeze()-np.min(PEA,axis=0).squeeze()
            SAL_mean[icluster,iexp]=np.mean(SAL,axis=0).squeeze()
            
#%%  
    vmin=0
    vmax=800                  
    iexp=0        
    Title='PEA annual max (Jm$^{-3}$)'
    Figname='../Figures/'+names[iexp]+'_PEA_max.png'    
    
    
    cluster_plot(x,y,PEA_max,vmin,vmax,Title,Figname,iexp=iexp)
    Title='PEA annual cycle (Jm$^{-3}$)'
    Figname='../Figures/'+names[iexp]+'_PEA_ann.png'
    cluster_plot(x,y,PEA_ann,vmin,vmax/2,Title,Figname,iexp=iexp)
#%%
    DPEA_max={}
    DPEA_ann={}
    DSAL_mean={}

    for icluster in range(clusters.values.shape[0]):
        DPEA_max[icluster]=PEA_max[icluster,0]-PEA_max[icluster,2]
        DPEA_ann[icluster]=PEA_ann[icluster,0]-PEA_ann[icluster,2]
        DSAL_mean[icluster]=SAL_mean[icluster,0]-SAL_mean[icluster,2]
        
#%%
    Title='DPEA annual max (Jm$^{-3}$)'
    vmin=-100
    vmax=100
    Figname='../Figures/'+names[0]+'_'+names[2]+'_DPEA_max.png'
    cluster_plot(x,y,DPEA_max,vmin,vmax,Title,Figname)
#%%
    Title='SAL annual mean'
    vmin=30
    vmax=36
    Figname='../Figures/'+names[0]+'_'+names[2]+'_SAL_mean.png'
    cluster_plot(x,y,SAL_mean,vmin,vmax,Title,Figname,iexp=0)

#%%
    Title='DSAL annual mean'
    vmin=-1
    vmax=1
    Figname='../Figures/'+names[0]+'_'+names[2]+'_DSAL_mean.png'
    cluster_plot(x,y,DSAL_mean,vmin,vmax,Title,Figname)
        