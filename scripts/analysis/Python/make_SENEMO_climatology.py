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
sys.path.insert(0,'/home/n01/n01/jholt/Git/COAsT/')
#needs branch     feature/535_stratification_diag
import coast


#Specify years to average
ystart=1990
ystop=2019

#ystop=2016
#ystop=1991
#EXPNAMS=['EXP_MES',  'EXP_MES_WAV',  'EXP_MES_WAV_NTM'] 
#EXPNAMS=[ 'EXP_MES_WAV_NTM_RIV']
names,dpaths,DOMS,_  = coast.experiments(experiments='experiments_paper.json')
for i,EXPNAM in enumerate(names):
    
    print(EXPNAM)
    print(names[i])
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
    #domain_outpath='/work/n01/n01/jholt/SENEMO/ASSESSMENT/'
    #nemo_out=coast.Gridded(fn_domain = fn_nemo_dom, config=fn_config_t_grid) #nemo_out = nemo_out. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))  
    DOMNAM='ORCA025-SE-NEMO'
    z_max=200
    fn_out='{0}/{1}/{1}_{2}_{3}_{4}_SST_SSS_PEA_MonClimate.nc'.format(domain_outpath,DOMNAM,ystart,ystop,EXPNAM)
        
    #Do the hardwork
    nemo_out=coast.GriddedMonthlyHydrographicClimatology(nemo,z_max=z_max)

    nemo_out.calc_climatologies()    
    #Write out as netcdf
    nemo_out.dataset.to_netcdf(fn_out)
