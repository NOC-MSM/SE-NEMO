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
def is_leap(year):
    """
    Checks if a given year is a leap year.

    Args:
        year: The year to check (integer).

    Returns:
        True if the year is a leap year, False otherwise.
    """
    # Rule 1: Must be divisible by 4
    if (year % 4) == 0:
        # Rule 2: If divisible by 100, must also be divisible by 400
        if (year % 100) == 0:
            if (year % 400) == 0:
                return True
            else:
                return False
        else:
            return True
    else:
        return False

names,dpaths,DOMS,_,year_start,year_stop  = coast.experiments(experiments='experiments_FC.json')

grid='T'
for i in [0,1,2,3]:#,1]:
    EXPNAM = names[i]
    ystart=year_start[i]
    ystop=year_stop[i]
    print(EXPNAM)


    domain_datapath=dpaths[i]   
    #make list of filenames
    fn_nemo_dat= coast.nemo_filename_maker(domain_datapath,ystart,ystop)            
    fn_nemo_dat=[]
    if 'CNRM' in EXPNAM:
        ESM='CNRM'
    else:
        ESM='GFDL'
    if 'hist' in EXPNAM:
        SSP='hist'
    else:
        SSP='ssp370'
    if 'bgc' in EXPNAM:
        variables = ['N3_n','P1_c','P2_c','P3_c','P4_c']
    else:
        variables = ['temperature', 'salinity']

    for year in range(ystart, ystop + 1):
        if 'bgc' in EXPNAM:
                new_name=f"{domain_datapath}/SE_{ESM}_subBGC_{SSP}_{year}.nc"
        else:
            days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            if is_leap(year):
                days[1]=29
            for Month in range(1, 13):
                day_stop=days[Month-1]
                new_name = f"{domain_datapath}/{EXPNAM}_1m_{year}{Month:02}01_{year}{Month:02}{day_stop:02}_grid_{grid}_{year}{Month:02}-{year}{Month:02}.nc"
                fn_nemo_dat.append(new_name)

    #Provide a domain.cfg file
    fn_nemo_dom=DOMS[i]
    
    #Provide a config file
    fn_config_t_grid='../Config/senemo_grid_t.json'    
                        
    #input datasets
    nemo = coast.Gridded(fn_data= fn_nemo_dat, fn_domain = fn_nemo_dom, config=fn_config_t_grid,multiple=True);#nemo = nemo.subset_as_copy(y_dim=range(860,1000),x_dim=range(1080,1180))
    #fix for nasty bug 
    nemo_dom=coast.Gridded(fn_domain = fn_nemo_dom, config=fn_config_t_grid);#nemo_dom = nemo_dom.subset_as_copy(y_dim=range(860,1000),x_dim=range(1080,1180))
    nemo.dataset['e3_0']=nemo_dom.dataset['e3_0']
    #Place to output data
    domain_outpath='/home/users/jholt/work/SENEMO/'

    DOMNAM='SENEMO_FC'
    z_max=200
    fn_out='{0}/{1}/{1}_NS_{2}_{3}_{4}_MonClimate.nc'.format(domain_outpath,DOMNAM,ystart,ystop,EXPNAM)
        
    #Do the hardwork
    nemo_out=coast.GriddedMonthlyHydrographicClimatology(nemo,z_max=z_max,variables=variables)

    nemo_out.calc_climatologies()    
    #Write out as netcdf
    nemo_out.dataset.to_netcdf(fn_out)
