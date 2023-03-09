#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 14:44:39 2022

@author: jholt
"""
# Needs SE-NEMO branch of coast
import socket
isliv = 'livljobs' in socket.gethostname()

import sys
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
import coast

#Specify years to average
ystart=1990
ystop=2019

#EXPNAMS=['EXP_MES',  'EXP_MES_WAV',  'EXP_MES_WAV_NTM'] 
#EXPNAMS=[ 'EXP_MES_WAV_NTM_RIV']
names,dpaths,DOMS,_  = coast.experiments(experiments='experiments_ZPS.json') 
for i,EXPNAM in enumerate(names):     
    print(EXPNAM)
#    domain_datapath='/gws/nopw/j04/class_vol2/senemo/jdha/SHORT_TESTS/' + EXPNAM +'/'
#    domain_datapath='/gws/nopw/j04/class_vol2/senemo/slwa/' + EXPNAM +'/'
    domain_datapath=dpaths[i]   
    #make list of filenames
    fn_nemo_dat= coast.nemo_filename_maker(domain_datapath,ystart,ystop)            
    
    #Provide a domain.cfg file
    fn_nemo_dom=DOMS[i]
    
    #Provide a config file
    fn_config_t_grid='../Config/senemo_grid_t.json'    
                        
    #input datasets
    nemo = coast.Gridded(fn_data= fn_nemo_dat, fn_domain = fn_nemo_dom, config=fn_config_t_grid,multiple=True);#nemo = nemo. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))
    #fix for nasty bug 
    nemo_dom=coast.Gridded(fn_domain = fn_nemo_dom, config=fn_config_t_grid) #;nemo_dom = nemo_dom. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))        
    nemo.dataset['e3_0']=nemo_dom.dataset['e3_0']
    #Place to output data
    domain_outpath='/home/users/jholt/work/SENEMO/ASSESSMENT/'
    
    nemo_out=coast.Gridded(fn_domain = fn_nemo_dom, config=fn_config_t_grid) #nemo_out = nemo_out. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))  
    DOMNAM='ORCA025-SE-NEMO'
    fn_out='{0}/{1}/{1}_{2}_{3}_{4}_SST_SSS_PEA_MonClimate.nc'.format(domain_outpath,DOMNAM,ystart,ystop,EXPNAM)
        
    #Do the hardwork
    coast.GriddedMonthlyHydrographicClimatology(nemo,nemo_out,Zmax=200)        
    #Write out as netcdf
    nemo_out.dataset.to_netcdf(fn_out)