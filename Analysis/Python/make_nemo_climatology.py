#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:27:35 2019

@author: jholt
"""
import socket
isliv = 'livljobs' in socket.gethostname()
import sys
import numpy as np
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/work/Git/COAsT/')    
import coast
Zmax=200
def calc_pea(nemo,Zd_mask):
#%%
 g=9.81
 Z=nemo.dataset.variables['depth_0'].values
 DZ=nemo.dataset.variables['e3_0'].values
 DP=np.sum(DZ*Zd_mask,axis=0)
 nz=nemo.dataset.dims['z_dim']
 nt=nemo.dataset.dims['t_dim']
 #if density and density_bar not in nemo, calculate them 
 if not 'density' in list(nemo.dataset.keys()):     
     nemo.construct_density(CT_AS=True,pot_dens=True)
 if not 'density_bar' in list(nemo.dataset.keys()):
     nemo.construct_density(CT_AS=True,rhobar=True,Zd_mask=Zd_mask,pot_dens=True)
 
 rho=nemo.dataset.variables['density'].values
 rho[np.isnan(rho)]=0
 rhobar=nemo.dataset.variables['density_bar'].values
 
 z_axis=0
 if len(nemo.dataset['density'].shape) == 4:   # includes time as first axis
  Z=np.repeat(Z[np.newaxis,:,:,:],nt,axis=0)
  DZ=np.repeat(DZ[np.newaxis,:,:,:],nt,axis=0)
  DP=np.repeat(DP[np.newaxis,:,:],nt,axis=0)  
  Zd_mask=np.repeat(Zd_mask[np.newaxis,:,:,:],nt,axis=0)  
  rhobar=np.repeat(rhobar[:,np.newaxis,:,:],nz,axis=1)  
  z_axis=1
 else:
  rhobar=np.repeat(rhobar[np.newaxis,:,:],nz,axis=0)  
 
 PEA=np.sum(Z*(rho-rhobar)*DZ*Zd_mask,axis=z_axis)*g/DP
#%%
 return PEA
#%%
def calc_zdmask(nemo,nemo_w,Zmax):
 """
    Calculates a 3D mask so a specified level Zmax. 1 for sea; 0 for below sea bed
    and linearly ramped for last level
 """
 Z=nemo.dataset.variables['depth_0'].values
 ZW=nemo_w.dataset.variables['depth_0'].values
 mbot=nemo.dataset.variables['bottom_level'].values
 mask=mbot!=0
 ZZ=ZW[1:,:,:]
 ZZ[ZZ==0]=np.nan
 ZW[1:,:,:]=ZZ

 Zd_mask=np.zeros((Z.shape))
 nz,ny,nx=np.shape(Z)
 kmax=np.zeros((ny,nx)).astype(int)
 IIkmax=np.zeros(np.shape(Z))
#
#careful assumes mbot is 1st sea point above bed ie new definition
 for i in range(nx):
     for j in range(ny):
       if mask[j,i]==1:          
         Zd_mask[0:mbot[j,i],j,i]=1 # mbot is not python style index so no +1
         kmax[j,i]=mbot[j,i]
         if ZW[mbot[j,i],j,i]>Zmax:
           kkmax=np.max(np.where(ZW[:,j,i]<Zmax))
           Zd_mask[kkmax+1:,j,i]=0
           Zd_mask[kkmax,j,i]=(Zmax-ZW[kkmax,j,i])/(ZW[kkmax+1,j,i]-ZW[kkmax,j,i])
           kmax[j,i]=kkmax
           IIkmax[kkmax,j,i]=1
 Ikmax=np.nonzero(IIkmax.ravel())          

 return Zd_mask,kmax,Ikmax
#%%%         

def make_climatology(nemo,nemo_w,DOMNAM,domain_outpath):
#%%
 Z=nemo.dataset.variables['depth_0'] 
 DZ=nemo.dataset.variables['e3_0']

 try:     
  A=np.load(domain_outpath + '/' +DOMNAM + '_Zd_mask.npz')
  
  Zd_mask=A['Zd_mask']
  kmax=A['kmax']
  Ikmax=A['kmax']
  print('read mask')
 except:
    print('calculating mask')
    Zd_mask,kmax,Ikmax=calc_zdmask(nemo,nemo_w,Zmax)
    np.savez(domain_outpath + '/' +DOMNAM + '_Zd_mask.npz',Zd_mask=Zd_mask,kmax=kmax,Ikmax=Ikmax)
#%%    
 ny=nemo.dataset.dims['y_dim']
 nx=nemo.dataset.dims['x_dim']
 nz=nemo.dataset.dims['z_dim']
 nt=nemo.dataset.dims['t_dim']
#%%
 print(nx,ny)
 SSTy=np.zeros((12,ny,nx))
 SSSy=np.zeros((12,ny,nx))
 PEAy=np.zeros((12,ny,nx))
 #NBTy=np.zeros((12,ny,nx))
 
#%%  
 tmp=nemo.dataset.variables['temperature']
 sal=nemo.dataset.variables['salinity']
 #PEA=calc_pea(nemo,Zd_mask)
 
 #need to find efficient method for bottom temperature
 #NBT=np.zeros((nt,ny,nx))
 #for it in range(nt): 
 #    NBT[it,:,:]=np.reshape(tmp[it,:,:,:].values.ravel()[Ikmax],(ny,nx))
     
 SST=tmp[:,0,:,:]
 SSS=sal[:,0,:,:]
 for im in range(12):
   print('Month',im)
   it=np.arange(im,nt,12).astype(int)
   SSTy[im,:,:]=np.mean(SST[it,:,:],axis=0)
   SSSy[im,:,:]=np.mean(SSS[it,:,:],axis=0)
#   PEAy[im,:,:]=np.mean(PEA[it,:,:],axis=0)
  # NBTy[im,:,:]=np.mean(NBT[it,:,:],axis=0)

 return SSTy,SSSy#,PEAy #,NBTy
#%%
def NEMO_FileNames(dpath,runtype,ystart,ystop):
#produce a list of nemo filenames
    names=[]    
    if runtype== 'SENEMO':
     for iy in range(ystart,ystop+1):
        for im in range(1,12+1):    
            MNTH=str(im);
            if im<10:
                 MNTH='0'+ MNTH
            YEAR=str(iy)
            new_name="{0}/SENEMO_1m_{1}0101_{1}1231_grid_T_{1}{2}-{1}{2}.nc".format(dpath,YEAR,MNTH)
            names.append(new_name)
    return names         
#%%    
if isliv: #NOC-liverpool
    domain_datapath='/work/jholt/JASMIN//SENEMO/NOTIDE/'
    domain_outpath='/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT/'
    domain_path='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/'
else: #JASMIN
    domain_datapath='/gws/nopw/j04/class_vol2/senemo/cwilso01/senemo/EXP_REF_NOTIDE/means/monthly/'
    domain_outpath='/home/users/jholt/work/SENEMO/'
    domain_path='/gws/nopw/j04/class_vol2/senemo/cwilso01/senemo/EXP_REF_NOTIDE/'

fn_nemo_dom=domain_path+'domcfg_eORCA025_v2.nc'
fn_config_t_grid='../Config/senemo_grid_t.json'
fn_nemo_dat=domain_datapath+'/SENEMO_1m_19800101_19801231_grid_T_1980*-1980*.nc'

#make list of filenames
ystart=1980
ystop=2009
ystop=1981
fn_nemo_dat= NEMO_FileNames(domain_datapath,'SENEMO',ystart,ystop)
  
fn_nemo_dom=domain_path+'domcfg_eORCA025_v2.nc'
fn_config_t_grid='/vkamino/work/jholt/Git/SE-NEMO/Analysis/Config/senemo_grid_t.json'    
#fn_nemo_dat=domain_datapath+'/SENEMO_1m_19800101_19801231_grid_T_1980*-1980*.nc'
#input datasets
nemo = coast.Gridded(fn_data= fn_nemo_dat, fn_domain = fn_nemo_dom, config=fn_config_t_grid,multiple=True)
nemo_w=coast.Gridded(fn_domain = fn_nemo_dom ,config='../Config/example_nemo_grid_w.json')
DOMNAM='ORCA025-SE-NEMO'
print('running')
#%%
SSTy,SSSy   = make_climatology(nemo,nemo_w,DOMNAM,domain_outpath)

 
