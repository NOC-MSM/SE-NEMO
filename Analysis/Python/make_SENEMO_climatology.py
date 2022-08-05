#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 14:44:39 2022

@author: jholt
"""
# Needs SE-NEMO branch of coast
import coast
#Specify years to average
ystart=1979
ystop=1979

EXPNAMS=['EXP_MES',  'EXP_MES_WAV',  'EXP_MES_WAV_NTM'] 
for EXPNAM in EXPNAMS:
    domain_datapath='/gws/nopw/j04/class_vol2/senemo/jdha/SHORT_TESTS/' + EXPNAM +'/'
    
    #make list of filenames
    fn_nemo_dat= coast.nemo_filenames(domain_datapath,'SENEMO',ystart,ystop)            
    
    #Provide a domain.cfg file
    fn_nemo_dom='/gws/nopw/j04/class_vol2/senemo/jdha/SHORT_TESTS/EXP_MES/domain_cfg_r015-r010_007_004v2.nc'
    
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
    coast.Annual_Climatology(nemo,nemo_out,Zmax=200)        
    #Write out as netcdf
    nemo_out.dataset.to_netcdf(fn_out)