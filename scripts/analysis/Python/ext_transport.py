#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 12:34:54 2023

@author: jholt
"""
import pickle
import socket
isliv = 'livljobs' in socket.gethostname()

import sys
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
sys.path.insert(0,'/home/n01/n01/jholt/Git/COAsT/')
import coast

import numpy as np
ystart=1976
ystop=2015
names,dpaths,DOMS,_  = coast.experiments(experiments='experiments_arch.json')
Q={}
for i,EXPNAM in enumerate(names):    

    print(EXPNAM)
#    domain_datapath='/gws/nopw/j04/class_vol2/senemo/jdha/SHORT_TESTS/' + EXPNAM +'/'
#    domain_datapath='/gws/nopw/j04/class_vol2/senemo/slwa/' + EXPNAM +'/'
    domain_datapath=dpaths[i]  
    fn_nemo_dat= coast.nemo_filename_maker(domain_datapath,ystart,ystop,grid='U')
    #Provide a config file
    fn_config_t_grid='../Config/senemo_grid_t.json'    
    imin=888
    imax=888+1
    jmin=330
    jmax=420                    
    #input datasets
    fn_nemo_dom=DOMS[i]
    nemo = coast.Gridded(fn_data= fn_nemo_dat, fn_domain = fn_nemo_dom, config=fn_config_t_grid,multiple=True);#nemo = nemo. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))
    nemo.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    #nemo.dataset.to_netcdf(f'{EXPNAM}_1174_275-428.nc')
    q=(nemo.dataset.uo * nemo.dataset.thkcello).sum(dim='depthu').values
    dy=nemo.dataset.e2.values
    DY=np.repeat(dy[np.newaxis,:,:],q.shape[0],axis=0)
    Q[EXPNAM]=np.sum(q*DY,axis=1).squeeze()
with open('so_flx','wb') as f:
    pickle.dump(Q, f)
    
