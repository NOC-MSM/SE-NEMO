#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 17:07:31 2022

@author: jholt
"""

import numpy as np
import matplotlib.pylab as plt
import scipy.io
import pandas as pd
import xarray as xr
import re
import sys
import socket
isliv = 'livljobs' in socket.gethostname()
if isliv:
    coast_dir='/login/jholt/work/Git/COAsT/'
else:
    coast_dir='/home/users/jholt/Git/COAsT/'
sys.path.insert(0,coast_dir)    
import coast
import matplotlib.pylab as plt
class CurrentsonT(coast.Gridded):
    """
    
    """
    def __init(self,fn_domain=None,config=None,**kwargs):
        gridded=coast.Gridded(fn_domain=fn_domain,config=config,**kwargs)
        self.dataset = gridded.dataset
    def currents_on_T(self,nemo_u,nemo_v):
#%%

        #U velocity on T points    
        UT=np.zeros(nemo_u.dataset.u_velocity.shape)            
        UT[:,:,:,1:]=0.5*(nemo_u.dataset.u_velocity.values[:,:,:,1:]
                         +nemo_u.dataset.u_velocity.values[:,:,:,:-1])
        #V velocity on T points
        VT=np.zeros(nemo_v.dataset.v_velocity.shape)
        VT[:,:,1:,:]=0.5*(nemo_v.dataset.v_velocity.values[:,:,1:,:]
                         +nemo_v.dataset.v_velocity.values[:,:,:-1,:])

        speed=np.sqrt(UT*UT+VT*VT)
        
        dims = nemo_u.dataset.u_velocity.dims        
        self.dataset["ut_velocity"]= xr.DataArray(UT,dims = dims)
        self.dataset["vt_velocity"]= xr.DataArray(VT,dims = dims )    
        self.dataset["speed_t"]= xr.DataArray(speed,dims = dims)  
#%%        
LME_Clusters='../Data/LME_Clusters_eORCA025.csv'
clusters = pd.read_csv(LME_Clusters)

for icluster in [1]: #range (6,17):
#%%    
    lims=np.array(clusters.values[icluster,2:6],dtype=int)
    datadir='/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/output/'
    domnam='/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/config/domain_cfg.nc'
    
    nemo_t=coast.CurrentsonT(fn_domain=domnam,config='example_nemo_grid_t.json',multiple =True,lims=lims)
    
    fn_data_u=[datadir+'SENEMO_1m_20070101_20071231_grid_U_200711-200711.nc']
    #           datadir+'SENEMO_1m_20070101_20071231_grid_U_200712-200712.nc'
               
    fn_data_v=[datadir+'SENEMO_1m_20070101_20071231_grid_V_200711-200711.nc']
    #           datadir+'SENEMO_1m_20070101_20071231_grid_V_200712-200712.nc'
               
    nemo_u=coast.Gridded(fn_domain=domnam,config='example_nemo_grid_u.json',
                         fn_data = fn_data_u, multiple =True,lims=lims)
    
    nemo_v=coast.Gridded(fn_domain=domnam,config='example_nemo_grid_v.json',
                         fn_data = fn_data_v, multiple =True,lims=lims)
    
    #nemo_t.subset(y_dim=range(lims[2],lims[3]),x_dim=range(lims[0],lims[1]))
    #nemo_u.subset(y_dim=range(lims[2],lims[3]),x_dim=range(lims[0],lims[1]))
    #nemo_v.subset(y_dim=range(lims[2],lims[3]),x_dim=range(lims[0],lims[1]))
    #%%
    nemo_t.subset(z_dim=[0])
    nemo_u.subset(z_dim=[0])
    nemo_v.subset(z_dim=[0])
    nemo_t.currents_on_T(nemo_u,nemo_v)
    plt.figure(icluster+1)
#    plt.pcolormesh(nemo_t.dataset.speed_t[0,0,:,:])
    name = clusters.values[icluster,1]
    nemo_t.plot_surface_circulation(name
                                  ,Vmax=0.5,Np=5
                                  ,headwidth=4,scale=50  )    
#nemo_t=currents_on_T(nemo_t,nemo_u,nemo_v)

