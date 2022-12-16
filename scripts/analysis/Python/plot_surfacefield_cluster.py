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
sys.path.insert(0, 'C:\\Users\\Jason.Goliath\\Documents\\GitHub\\COAsT\\')
import coast
import matplotlib.pylab as plt

LME_Clusters='../Data/LME_Clusters_eORCA025.csv'
clusters = pd.read_csv(LME_Clusters)


LME_Data=np.load('../Data/LME_gridinfo_V4.npz')
a=scipy.io.loadmat('../Data/ORCA025_ROAM_GLB_LMEmaskV4.mat')
J_offset=186
LME_mask=a['LME_mask'][:,:].T

bathyname='../Data/eORCA025_bathy_meter.nc'
config='example_nemo_grid_t.json'
bathy=coast.Gridded(fn_data=bathyname,config=config)
lon=bathy.dataset.longitude
lat=bathy.dataset.latitude


icluster=15
lims=clusters.values[icluster,2:6]

lims=np.array(clusters.values[icluster,2:6],dtype=int)
Lims=np.copy(lims)
Lims[2]=lims[2]-J_offset
Lims[3]=lims[3]-J_offset

M=LME_mask[Lims[2]:Lims[3]+1,Lims[0]:Lims[1]+1]
#%%
I=np.where(~pd.isnull(clusters.values[icluster,6:]))[0]
LMEs=[]
for i in I:
    LMEs.append(int(re.findall(r'\d+',clusters.values[icluster,6+i])[0]))
MM=np.ones(M.shape)*np.nan
for LME in LMEs:
    MM[M==LME]=LME

Bathy=bathy.subset_as_copy(y_dim=range(lims[2],lims[3]),x_dim=range(lims[0],lims[1]))
D=np.ma.masked_where(Bathy.dataset.Bathymetry.values==0,Bathy.dataset.Bathymetry.values)

    
#%%




plt.pcolormesh(MM)