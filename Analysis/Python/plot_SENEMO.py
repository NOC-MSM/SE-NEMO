#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 09:42:55 2021

@author: jholt
"""
import sys
sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
import matplotlib.pylab as plt
import coast
fn_nemo_dom='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc'
fn_nemo_dat1='/work/jholt/JASMIN//SENEMO/NOTIDE/SENEMO_1m_20000101_20001231_grid_T_200012-200012.nc'
fn_nemo_dat2=  '/work/jholt/JASMIN//SENEMO/TIDE/SENEMO_1m_20000101_20001231_grid_T_200012-200012.nc'

fn_config_t_grid='/login/jholt/work/Git/COAsT/config/example_nemo_grid_t.json'
nemo_t1 = coast.Gridded(fn_data = fn_nemo_dat1, fn_domain = fn_nemo_dom, config=fn_config_t_grid)
nemo_t2 = coast.Gridded(fn_data = fn_nemo_dat2, fn_domain = fn_nemo_dom, config=fn_config_t_grid)

SSS1=nemo_t1.dataset.variables['so_abs'].values[0,0,:,:]
SSS2=nemo_t2.dataset.variables['so_abs'].values[0,0,:,:]
SST1=nemo_t1.dataset.variables['thetao_con'].values[0,0,:,:]
SST2=nemo_t2.dataset.variables['thetao_con'].values[0,0,:,:]


plt.pcolormesh(SSS2)