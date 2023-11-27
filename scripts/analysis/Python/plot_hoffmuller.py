#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 12:34:54 2023

@author: jholt
"""

import socket
isliv = 'livljobs' in socket.gethostname()

import sys
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
sys.path.insert(0,'/home/n01/n01/jholt/Git/COAsT/')
import coast
import matplotlib.pylab as plt
ystart=1977
ystop=2011
names,dpaths,DOMS,_  = coast.experiments(experiments='experiments_arch.json')
for i,EXPNAM in enumerate([names[1]]):     
    print(EXPNAM)
#    domain_datapath='/gws/nopw/j04/class_vol2/senemo/jdha/SHORT_TESTS/' + EXPNAM +'/'
#    domain_datapath='/gws/nopw/j04/class_vol2/senemo/slwa/' + EXPNAM +'/'
    domain_datapath=dpaths[i]  
    fn_nemo_dat= coast.nemo_filename_maker(domain_datapath,ystart,ystop)
    #Provide a config file
    fn_config_t_grid='../Config/senemo_grid_t.json'    
    imin=1174
    imax=1174+1
    jmin=275
    jmax=428                    
    #input datasets
    fn_nemo_dom=DOMS[i]
    nemo = coast.Gridded(fn_data= fn_nemo_dat, fn_domain = fn_nemo_dom, config=fn_config_t_grid,multiple=True);#nemo = nemo. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))
    nemo.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax),z_dim=[0])