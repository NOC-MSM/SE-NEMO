#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:47:15 2023

@author: jholt
"""

import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
import coast
import surfacefields as sf
cmap1=sf.lightcolormap(32,2)
cmap1.set_bad([0.75,0.75,0.75])
x_min=-28;x_max=12;y_min=40;y_max=69
x_min=-98;x_max=26.5;y_min=-56;y_max=69
config='example_nemo_grid_t.json'
names,dpaths,DOMS,_  = coast.experiments(experiments='experiments.json')
PEA_max={}
for iexp in [0,2]:
    EXP_NAM=names[iexp]
    
    fn_data='/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/ORCA025-SE-NEMO_1990_2019_'+EXP_NAM+'_SST_SSS_PEA_MonClimate.nc'
    
    nemo_t=coast.Gridded(fn_data=fn_data,config=config)

    j,i,_=nemo_t.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
    imin=min(i)
    imax=max(i)
    jmin=min(j)
    jmax=max(j) 
    nemo_t.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    PEA_max[iexp]=np.ma.masked_where(nemo_t.dataset.bottom_level.values==0,
                                     np.max(nemo_t.dataset.PEA_monthy_clim.values,axis=0)).squeeze()


#%%
fig,axs=plt.subplots(nrows=1,ncols=2,figsize=[11.69,8.27])

im=axs[0].pcolormesh(PEA_max[0],vmin=0,vmax=800,cmap=cmap1)
axs[0].set_xticks([]);axs[0].set_yticks([])
axs[0].set_title('GS1p2_full PEA annual max 1990-2019 (Jm$^{-3}$)')
fig.colorbar(im,ax=axs[0],orientation='vertical')
im2=axs[1].pcolormesh(PEA_max[0]-PEA_max[2],vmin=-100,vmax=100,cmap=cmap1)
axs[1].set_xticks([]);axs[1].set_yticks([])
fig.colorbar(im2,ax=axs[1],orientation='vertical')
axs[1].set_title('GS1p2 - GS1p0 PEA annual max  (Jm$^{-3}$)')


plt.savefig('../Figures/PEAmax_Atlantic_GS1p2_GS1p1.png',dpi=300)