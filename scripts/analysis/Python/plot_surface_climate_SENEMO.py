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
plt.ion()
cmap1=sf.lightcolormap(32,2)
cmap1.set_bad([0.75,0.75,0.75])
x_min=-28;x_max=12;y_min=40;y_max=69
x_min=-98;x_max=26.5;y_min=-56;y_max=69
config='example_nemo_grid_t.json'

names,dpaths,DOMS,_,year_start,year_stop  = coast.experiments(experiments='experiments_FC.json')
netppmean={}
for iexp in [2,3]:
    EXPNAM=names[iexp]
    ystart=year_start[iexp]
    ystop=year_stop[iexp]
    domain_outpath='/home/users/jholt/work/SENEMO/'

    DOMNAM='SENEMO_FC'
    z_max=200
    fn_data='{0}/{1}/{1}_NS_{2}_{3}_{4}_MonClimate.nc'.format(domain_outpath,DOMNAM,ystart,ystop,EXPNAM)
    fn_domain=DOMS[iexp]
    nemo_t=coast.Gridded(fn_data=fn_data,fn_domain=fn_domain,config=config, no_depths=True)

    j,i,_=nemo_t.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
    imin=min(i)
    imax=max(i)
    jmin=min(j)
    jmax=max(j) 
    nemo_t.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    netppmean[iexp]=np.ma.masked_where(nemo_t.dataset.bottom_level.values==0,
                                     np.mean(nemo_t.dataset['Ptot_NPP_result_monthly_clim'].values,axis=0)
                                       *365/1000
                                       ).squeeze()


#%%
axs=[0,0]
fig,axs[1]=plt.subplots(nrows=1,ncols=1,figsize=[11.69,8.27])

#im=axs[0].pcolormesh(netppmean[2],vmin=0,vmax=200,cmap=cmap1)
#axs[0].set_xticks([]);axs[0].set_yticks([])
#axs[0].set_title('netPP surface mean 1991-2000 gCm$^{3}$yr$^{-1}$')
#fig.colorbar(im,ax=axs[0],orientation='vertical')
im2=axs[1].pcolormesh(netppmean[3]/netppmean[2]-1,vmin=-1,vmax=1,cmap=cmap1)
axs[1].set_xticks([]);axs[1].set_yticks([])
fig.colorbar(im2,ax=axs[1],orientation='vertical')
axs[1].set_title('$\Delta$netPP surface mean (2061-2070) - (1991-2000) Fraction ')


plt.savefig(f'../Figures/netPP_{EXPNAM}.png',dpi=300)