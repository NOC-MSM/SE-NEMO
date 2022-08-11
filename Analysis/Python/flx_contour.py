#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 16:46:08 2022

@author: jholt
"""

import sys
sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
import matplotlib.pylab as plt
import coast
import numpy as np
fn_config_t_grid='/login/jholt/work/Git/COAsT/config/example_nemo_grid_t.json'
fn_config_f_grid='/login/jholt/work/Git/COAsT/config/example_nemo_grid_f.json'
fn_config_u_grid='/login/jholt/work/Git/COAsT/config/example_nemo_grid_u.json'
fn_config_v_grid='/login/jholt/work/Git/COAsT/config/example_nemo_grid_v.json'
    
DOMS=['/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc',
      '/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc',
      '/work/jholt/JASMIN//SENEMO/JDHA/SHORT_TESTS/EXP_MES_WAV_NTM/domain_cfg_r015-r010_007_004v2.nc',
      '/work/jholt/JASMIN//SENEMO/JDHA/TO_SORT/SE-NEMO/EXP_SZT39_TAPER_TIDE/domain_cfg_ztaper_match.nc'
      ]
#DOMS=['/work/jholt/JASMIN//SENEMO/JDHA/TO_SORT/SE-NEMO/EXP_SZT39_TAPER_TIDE/domain_cfg_ztaper_match.nc']
DOMS=['/work/jholt/JASMIN//SENEMO/JDHA/SHORT_TESTS/EXP_MES_WAV_NTM/domain_cfg_r015-r010_007_004v2.nc']      
#DOMS=['/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc']      

dpaths=['/work/jholt/JASMIN//SENEMO/NOTIDE/',
       '/work/jholt/JASMIN/SENEMO/TIDE/outputs/',
       '/work/jholt/JASMIN//SENEMO/JDHA/SHORT_TESTS/EXP_MES_WAV_NTM/',
       '/work/jholt/JASMIN//SENEMO/JDHA/LONG_TESTS/EXP_SZT_TIDE/'
       ]
#dpaths=['/work/jholt/JASMIN//SENEMO/JDHA/LONG_TESTS/EXP_SZT_TIDE/']
dpaths=['/work/jholt/JASMIN//SENEMO/JDHA/SHORT_TESTS/EXP_MES_WAV_NTM/']
#names=['NOTIDE','Tide','EXP_MES_WAV_NTM','EXP_SZT_TIDE']
names=['EXP_MES_WAV_NTM']
#names=['EXP_SZT_TIDE']
Unmean={}
Untm={}

Zmean={}
Un={}
Z={}
ystart=1981
ystop=1981
for i,name in enumerate(names): 

    fn_nemo_dat_u= coast.nemo_filenames(dpaths[i],'SENEMO',ystart,ystop,grid='U') 
    fn_nemo_dat_v= coast.nemo_filenames(dpaths[i],'SENEMO',ystart,ystop,grid='V') 
       
    fn_nemo_dom=DOMS[i]
#    fn_nemo_dat_u=dpaths[i]+'SENEMO_1m_19800101_19801231_grid_U_198007-198007.nc'
#    fn_nemo_dat_v=dpaths[i]+'SENEMO_1m_19800101_19801231_grid_V_198007-198007.nc'

    
    
    
    nemo_t = coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_t_grid,calc_bathy=True)    
    nemo_f = coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_f_grid,calc_bathy=True)
    nemo_u = coast.Gridded(fn_data=fn_nemo_dat_u, fn_domain=fn_nemo_dom, config=fn_config_u_grid,multiple=True)
    nemo_v = coast.Gridded(fn_data=fn_nemo_dat_v, fn_domain=fn_nemo_dom, config=fn_config_v_grid,multiple=True)
    
    nemo_f.subset(y_dim=range(860,1000),x_dim=range(1080,1180))
    nemo_u.subset(y_dim=range(860,1000),x_dim=range(1080,1180))
    nemo_v.subset(y_dim=range(860,1000),x_dim=range(1080,1180))
    
    contours, no_contours = coast.Contour.get_contours(nemo_f, 200)
    y_ind, x_ind, contour = coast.Contour.get_contour_segment(nemo_f, contours[0], [50, -10], [60, 2])
    cont_f = coast.ContourF(nemo_f, y_ind, x_ind, 200)
    cont_f.calc_cross_contour_flow(nemo_u, nemo_v)

    Un[name]=cont_f.data_cross_flow.normal_velocities.values
    Un[name][Un[name]==0]=np.nan
    Untm[name]=np.nanmean(Un[name],axis=0)
    Unmean[name]=np.nanmean(Untm[name],axis=1)
    Z[name]=-cont_f.data_cross_flow.depth_0.values
    Zmean[name]=np.nanmean(Z[name],axis=1)
    r_dim=cont_f.data_cross_flow.r_dim.values
#%%
plt.plot(Unmean[names[0]],Zmean[names[0]],'o-',
         Unmean[names[1]],Zmean[names[1]],'o-',
         Unmean[names[2]],Zmean[names[2]],'o-',
         Unmean[names[3]],Zmean[names[3]],'o-'         
         )
plt.legend(names)

