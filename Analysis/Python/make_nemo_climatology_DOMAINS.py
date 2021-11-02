#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 16:12:00 2019

@author: jholt
"""
import sys
sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
import numpy as np
import coast
import matplotlib.pylab as plt
import make_nemo_climatology as mc
import pickle as pkl
domain_datapath='/work/jholt/JASMIN//SENEMO/NOTIDE/'
domain_outpath='/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT/'
domain_path='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/'

domain_datapath='/gws/nopw/j04/class_vol2/senemo/cwilso01/senemo/EXP_REF_NOTIDE/means/monthly/'
domain_outpath='/home/users/jholt/work/SENEMO/'
domain_path='/gws/nopw/j04/class_vol2/senemo/cwilso01/senemo/EXP_REF_NOTIDE/'
DOMNAM='ORCA025-SE-NEMO'

fn_nemo_dom=domain_path+'domcfg_eORCA025_v2.nc'
fn_config_t_grid='/vkamino/work/jholt/Git/SE-NEMO/Analysis/Config/senemo_grid_t.json'
fn_nemo_dat=domain_datapath+'/SENEMO_1m_19800101_19801231_grid_T_1980*-1980*.nc'

DOMNAM='ORCA025-SE-NEMO'
  
fn_nemo_dom=domain_path+'domcfg_eORCA025_v2.nc'
fn_config_t_grid='/vkamino/work/jholt/Git/SE-NEMO/Analysis/Config/senemo_grid_t.json'    
fn_nemo_dat=domain_datapath+'/SENEMO_1m_19800101_19801231_grid_T_1980*-1980*.nc'
nemo = coast.Gridded(fn_data= fn_nemo_dat, fn_domain = fn_nemo_dom, config=fn_config_t_grid,multiple=True)
nemo_w=coast.Gridded(fn_domain = fn_nemo_dom ,config='../Config/example_nemo_grid_w.json')

print('running')
#%%

SSTy,SSSy, PEAy, NBTy =mc.make_climatology(nemo,nemo_w,DOMNAM,domain_outpath)

Mask=np.repeat(grd.mask[:,:,np.newaxis],12,axis=2)

outname= e.assess_path + '/' + DOMNAM  + '/' +DOMNAM + '_' + EXPNUM +'_TSclim_NOTIDE' + str(yearstart) + '_' +str(yearstop)+'.nc'
f = Dataset(outname, 'w',  format='NETCDF4')
f.createDimension('lon',grd.nx)
f.createDimension('lat',grd.ny)
f.createDimension('mnth',12)
Lon_dom=f.createVariable('lon_dom', np.float64, ('lon','lat'))
Lat_dom=f.createVariable('lat_dom', np.float64, ('lon','lat'))
ssty=f.createVariable('SSTy', np.float64, ('lon','lat','mnth'))
sssy=f.createVariable('SSSy', np.float64, ('lon','lat','mnth'))
peay=f.createVariable('PEAy', np.float64, ('lon','lat','mnth'))
nbty=f.createVariable('NBTy', np.float64, ('lon','lat','mnth'))

#SSTy=np.ma.masked_where(Mask,SSTy)

Lon_dom[:,:]=grd.lon_dom
Lat_dom[:,:]=grd.lat_dom
ssty[:,:,:]=SSTy
sssy[:,:,:]=SSSy
peay[:,:,:]=PEAy
nbty[:,:,:]=NBTy
f.close()
 
