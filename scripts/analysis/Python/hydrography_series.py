#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 15:23:41 2023

@author: jholt
"""

import socket
isliv = 'livljobs' in socket.gethostname()
import numpy as np
import matplotlib.pylab as plt
import sys
import pandas as pd

if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
 
import coast
import scipy.io

J_offset=186 #account for extra rows in eORCA if data is made for normal ORCA

A=np.load('../Data/LME_gridinfo_V4.npz')
a=scipy.io.loadmat('../Data/ORCA025_ROAM_GLB_LMEmaskV4.mat')
LME_mask=a['LME_mask'][:,:].T
if isliv:
 datapath_LME ='/work/jholt/Git/DEV_jholt/GCO_Hydro/Data/'
 fn_bathymetry='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/eORCA025_bathy_meter.nc'
else:
 datapath_LME ='/home/users/jholt/work/SENEMO/ASSESSMENT/EN4/'
 fn_bathymetry='/home/users/jholt/work/SENEMO/senemo/INPUTS/eORCA025_bathy_meter.nc'    
 Assessdir='/home/users/jholt//work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/'

DATANAME='ORCA025'
vnames=['PEA_monthy_clim', 'SST_monthy_clim','SSS_monthy_clim']
RUNNAMS=[
'ORCA025-SE-NEMO_1990_2019_EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2',
'ORCA025-SE-NEMO_1990_2019_ZPS_TIDE',
'ORCA025-SE-NEMO_1990_2019_ZPS_NOTIDE'  
    ]
PEAym=np.zeros((12,len(RUNNAMS)))
for iR,RUNNAM in enumerate(RUNNAMS):
    fn_nemo_dat=Assessdir +RUNNAM+'_SST_SSS_PEA_MonClimate.nc'
    f = coast.Gridded(fn_data= fn_nemo_dat)
    f_bathy=coast.Gridded(fn_data= fn_bathymetry)
    
    
    for ilme in [21]:
        LMENAM=A['DOMNAM'][ilme]
        i_min1=1134
        j_min1=960
        i_max1=1148
        j_max1=974
    

        Depth=f_bathy.dataset.variables['Bathymetry'].values[j_min1:j_max1+1,i_min1:i_max1+1]     
        PEAy=f.dataset.variables[vnames[0]].values[:,j_min1:j_max1+1,i_min1:i_max1+1]
        SSTy=f.dataset.variables[vnames[1]].values[:,j_min1:j_max1+1,i_min1:i_max1+1]
        SSSy=f.dataset.variables[vnames[2]].values[:,j_min1:j_max1+1,i_min1:i_max1+1]
        mask=f.dataset.variables['bottom_level'].values[j_min1:j_max1+1,i_min1:i_max1+1] !=0
        nx,ny=mask.shape
        PEAym[:,iR]=np.nanmean(PEAy.reshape(12,ny*nx),axis=1)
        EN4=coast.Gridded(datapath_LME + '/' + DATANAME  +'/'+ LMENAM + '_EN4mnthgrid_V3.nc')
        lims=EN4.dataset.lims.values
        i_min=i_min1-lims[0]
        i_max=i_min+i_max1-i_min1
        j_min=j_min1-lims[2]-J_offset
        j_max=j_min+j_max1-j_min1
        SST_EN4_LME=EN4.dataset.SST_g.values.T[:,j_min:j_max+1,i_min:i_max+1]
        SSS_EN4_LME=EN4.dataset.SSS_g.values.T[:,j_min:j_max+1,i_min:i_max+1]
        PEA_EN4_LME=EN4.dataset.PEA_g.values.T[:,j_min:j_max+1,i_min:i_max+1]

        PEA_EN4_ds=PEA_EN4_LME.reshape(12,ny*nx)
       
        pea_std=np.repeat(np.nanstd(PEA_EN4_ds,axis=1)[:,np.newaxis],nx*ny,axis=1)
        pea_mean=np.repeat(np.nanmean(PEA_EN4_ds,axis=1)[:,np.newaxis],nx*ny,axis=1)        

        PEA_EN4_ds[np.abs(PEA_EN4_ds-pea_mean)>2.*pea_std]=np.nan
        PEAm_EN4=np.nanmean(PEA_EN4_ds,axis=1)
t=range(1,13)
plt.plot(t,PEAm_EN4,'-o',t,PEAym)


