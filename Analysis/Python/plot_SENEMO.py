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
import numpy as np
fn_nemo_dom='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc'
fn_nemo_dat1='/work/jholt/JASMIN//SENEMO/NOTIDE/SENEMO_1m_19800101_19801231_grid_T_198001-198001.nc'
fn_nemo_dat2=  '/work/jholt/JASMIN//SENEMO/TIDE/SENEMO_1m_19800101_19801231_grid_T_198001-198001.nc'
fn_nemo_dat_w1='/work/jholt/JASMIN//SENEMO/NOTIDE/SENEMO_1m_19800101_19801231_grid_W_198001-198001.nc'
fn_nemo_dat_w2=  '/work/jholt/JASMIN//SENEMO/TIDE/SENEMO_1m_19800101_19801231_grid_W_198001-198001.nc'

fn_config_t_grid='/login/jholt/work/Git/COAsT/config/example_nemo_grid_t.json'
fn_config_w_grid='/login/jholt/work/Git/COAsT/config/example_nemo_grid_w.json'


nemo_t1 = coast.Gridded(fn_data = fn_nemo_dat1, fn_domain = fn_nemo_dom, config=fn_config_t_grid)
nemo_t2 = coast.Gridded(fn_data = fn_nemo_dat2, fn_domain = fn_nemo_dom, config=fn_config_t_grid)


nemo_w1 = coast.Gridded(fn_data = fn_nemo_dat_w1, fn_domain = fn_nemo_dom, config=fn_config_w_grid)
nemo_w2 = coast.Gridded(fn_data = fn_nemo_dat_w2, fn_domain = fn_nemo_dom, config=fn_config_w_grid)



SSS1=nemo_t1.dataset.variables['so_abs'].values[0,0,:,:]
SSS2=nemo_t2.dataset.variables['so_abs'].values[0,0,:,:]
SST1=nemo_t1.dataset.variables['thetao_con'].values[0,0,:,:]
SST2=nemo_t2.dataset.variables['thetao_con'].values[0,0,:,:]


plt.pcolormesh(SSS2)

#%%
j=660
i=200


j=950
i=1140

i=565
j=685

avh1=nemo_w1.dataset.variables['difvho'].values[0,:,j,i]
avh2=nemo_w2.dataset.variables['difvho'].values[0,:,j,i]
T1=nemo_t1.dataset.variables['thetao_con'].values[0,:,j,i]
T2=nemo_t2.dataset.variables['thetao_con'].values[0,:,j,i]

Z=nemo_t2.dataset.coords['depth_0'].values[:,j,i]



plt.subplot(1,2,1)
plt.plot(T1,-Z,T2,-Z,'--')
plt.legend(['No Tide tmp','Tide tmp'])
plt.title('Profiles at j,i={0},{1}'.format(j,i))
plt.subplot(1,2,2)
plt.plot(np.log10(avh1),-Z,np.log10(avh2),-Z,'--')
plt.legend(['No Tide avt (log10)','Tide avt (log10)'])

plt.savefig('/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT//ORCA025-SE-NEMO/prof_j{0}_i{1}_1980_tide-notide.png'.format(j,i))






























