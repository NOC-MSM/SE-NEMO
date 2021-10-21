#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:27:35 2019

@author: jholt
"""
import numpy as np
from netCDF4 import Dataset

import nemogrid
import read_nemo_field as rf
import environment as e
from  seawater import eos80 as sw
#Depth mask   
#Zd_mask[ip,0:int(kbot[ip])]=1
#Imax=np.max(np.nonzero(ZW>-Zmax)[0]) #note ZW index
#if Imax<kbot[ip]:
# Zd_mask[ip,Imax:nz]=0
# Zd_mask[ip,Imax]=(ZW[Imax]-(-Zmax))/(ZW[Imax]-ZW[Imax+1])
#if Zd_mask[ip,Imax-1] < 0 or Zd_mask[ip,Imax-1]  >1:
# print('error', ip,Zd_mask[ip,Imax-1]),Imax,kbot[ip]
#domain_outpath='/projectsa/Recicle/outputs/v2/'
#DOMNAM='HadGEM2-ES'
#EXPNUM=''
#RUNNAM=''
#
#yearstart=2024
#yearstop=2052
Zmax=200
def calc_pea(tmp,sal,Z, DZ,Zd_mask):
#not finiished yet    
 g=9.81
 nz=tmp.shape[2]
 DP=np.sum(DZ*Zd_mask,axis=2)
 Tbar=np.sum(tmp*DZ*Zd_mask,axis=2)/DP
 Sbar=np.sum(sal*DZ*Zd_mask,axis=2)/DP
 rho=sw.dens0(sal,tmp)
 rhobar=sw.dens0(Sbar,Tbar)
 rhobar_2d=np.repeat(rhobar[:,:,np.newaxis],nz,axis=2)
 PEA=np.sum(Z*(rho-rhobar_2d)*DZ*Zd_mask,axis=2)*g/DP
 return PEA

def calc_zdmask(nemo,nemo_w,Zmax):
#%%    
 Z=nemo.dataset.variables['depth_0']
 ZW=nemo_w.dataset.variables['depth_0']
 ZZ=ZW[:,:,1:]
 ZZ[ZZ==0]=np.nan
 ZW[:,:,1:]=ZZ
 mbot=np.squeeze(grd.mbot[:,:])
 Zd_mask=np.zeros((Z.shape))
 nx,ny,nz=np.shape(Z)
 kmax=np.zeros((nx,ny)).astype(int)
 IIkmax=np.zeros(np.shape(Z))
#%% 
#careful assumes mboit is 1st sea point above bed ie new definition
 for i in range(nx):
     for j in range(ny):
       if grd.mask[i,j]==1:          
         Zd_mask[i,j,0:mbot[i,j]]=1 # mbot is not python style index so no +1
         kmax[i,j]=mbot[i,j]
         if ZW[i,j,mbot[i,j]]>Zmax:
           kkmax=np.max(np.where(ZW[i,j,:]<Zmax))
           Zd_mask[i,j,kkmax+1:]=0
           Zd_mask[i,j,kkmax]=(Zmax-ZW[i,j,kkmax])/(ZW[i,j,kkmax+1]-ZW[i,j,kkmax])
           kmax[i,j]=kkmax
           IIkmax[i,j,kkmax]=1
 Ikmax=np.nonzero(IIkmax.ravel())          
#%%
 return Zd_mask,kmax,Ikmax
         

def make_climatology(nemo,nemo_w,yearstart,yearstop):
 Z=nemo.dataset.variables['depth_0'] 
 DZ=nemo.dataset.variables['e3_0']
 print('calculating mask')
 try:
#%%     
  A=np.load(e.assess_path + '/' + DOMNAM  + '/' +DOMNAM + '_Zd_mask.npz')
  Zd_mask=A['Zd_mask']
  kmax=A['kmax']
  Ikmax=A['kmax']
#%%  
 except:
    Zd_mask,kmax,Ikmax=calc_zdmask(nemo,nemo_w,Zmax)
    np.savez(e.assess_path + '/' + DOMNAM  + '/' +DOMNAM + '_Zd_mask.npz',Zd_mask=Zd_mask,kmax=kmax,Ikmax=Ikmax)
    
 nx,ny,nz=np.shape(Z)
 def make_climatology_monthly(nemo,yearstart,yearstop):
#%%
  nyear = yearstop-yearstart+1;
  print(nx,ny)
  SSTy=np.zeros((nx,ny,12))
  SSSy=np.zeros((nx,ny,12))
  PEAy=np.zeros((nx,ny,12))
  NBTy=np.zeros((nx,ny,12))
  #PEATy=np.zeros((nx,ny,12))
#%%  
  tmp=nemo.dataset.variables['temperature']
  sal=nemo.dataset.variables['salinity']
  rho=nemo.construct_density
  for iy in np.arange(yearstart,yearstop+1):    
   print(iy)
   for im in np.arange(0,12):
#%%
    print(iy,im,DOMNAM, domain_outpath,EXPNUM,RUNNAM )
#    tmp,sal =rf.read_nemo(iy,np.int(im)+1,domain_outpath,DOMNAM,EXPNUM,RUNNAM,'1m','1M','sosstsst','sosaline')       
    if FT=='SSB':
     FRQ='1d'
    else:
     FRQ='1m'   
    tmp,sal =rf.read_nemo(iy,np.int(im)+1,domain_outpath,DOMNAM,EXPNUM,RUNNAM,FRQ,FT,tmpvar,salvar)       
    print('read in')
    if FT=='SSB':
     tmp=np.mean(tmp,axis=3)
     sal=np.mean(sal,axis=3)
    else:    
     tmp=np.squeeze(tmp)
     sal=np.squeeze(sal)
    pea=calc_pea(tmp,sal,Z, DZ,Zd_mask)
    print('pea calc')
    SSTy[:,:,im]=SSTy[:,:,im]+tmp[:,:,0]
    SSSy[:,:,im]=SSSy[:,:,im]+sal[:,:,0]
    PEAy[:,:,im]=PEAy[:,:,im]+pea[:,:]
#this shoudl be faster:    
#    tmpr=tmp.ravel()[Ikmax]
#    NBTy[i,j,im]=NBTy[i,j,im]+tmpr
    for i in range(nx):
     for j in range(ny):   
      NBTy[i,j,im]=NBTy[i,j,im]+tmp[i,j,kmax[i,j]]
    print('NBT calc')
#    PEATy[:,:,im]=PEATy[:,:,im]+peat[:,:,im]
#%%    
  SSTy=SSTy/nyear
  SSSy=SSSy/nyear
  PEAy=PEAy/nyear
  NBTy=NBTy/nyear
  #PEATy=PEATy/nyear
#%%  
  return SSTy,SSSy,PEAy,NBTy


#tmpmr=tmpm.reshape(nx*ny,nz)
#salmr=salm.reshape(nx*ny,nz)

# make_climateology_daily_month(yearstart,yearstop,domain_outpath,DOMNAM,EXPNUM,RUNNAM)
  
 if FT=='M':
  SSTy,SSSy=make_climatology_daily_month(yearstart,yearstop,domain_outpath,DOMNAM,EXPNUM,RUNNAM)
  return SSTy,SSSy 
 elif  FT=='1M' or FT=='CMEMS' or FT=='SSB' or FT=='SENEMO':
  SSTy,SSSy, PEAy, NBTy =make_climatology_monthly(yearstart,yearstop,domain_outpath,DOMNAM,EXPNUM,RUNNAM,FT)  
  return SSTy,SSSy, PEAy, NBTy
 else:
  PEAy,SSTy,SSSy, PEATy=make_climatology_monthly_year(yearstart,yearstop,domain_outpath,DOMNAM,EXPNUM,RUNNAM)
  return PEAy,SSTy,SSSy , PEATy

if __name__=='__main__':
 domain_outpath='/projectsa/accord/GCOMS1k/OUTPUTS/'
 DOMNAM='GTHI35'
 EXPNUM='01'
 RUNNAM=''
 yearstart=1995
 yearstop=2015
 yearstop=2011
 grd=nemogrid.nemo_grid(DOMNAM,DOMNAM)
 SSTy,SSSy=make_climatology(yearstart,yearstop,domain_outpath,DOMNAM,EXPNUM,RUNNAM,'M')
 outname= e.assess_path + '/' + DOMNAM  + '/' +DOMNAM + '_' + EXPNUM +'_TSclim_' + str(yearstart) + '_' +str(yearstop)+'.nc'
 f = Dataset(outname, 'w',  format='NETCDF4')
 f.createDimension('lon',grd.nx)
 f.createDimension('lat',grd.ny)
 f.createDimension('mnth',12)
 Lon_dom=f.createVariable('lon_dom', np.float64, ('lon','lat'))
 Lat_dom=f.createVariable('lat_dom', np.float64, ('lon','lat'))
 ssty=f.createVariable('SSTy', np.float64, ('lon','lat','mnth'))
 sssy=f.createVariable('SSSy', np.float64, ('lon','lat','mnth'))
 Lon_dom[:,:]=grd.lon_dom
 Lat_dom[:,:]=grd.lat_dom
 ssty[:,:,:]=SSTy
 sssy[:,:,:]=SSSy
 f.close()
 
 
 
 

