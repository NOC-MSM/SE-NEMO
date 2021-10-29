#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:27:35 2019

@author: jholt
"""
import numpy as np

Zmax=200
def calc_pea(nemo,Zd_mask):
#%%
 g=9.81
 Z=nemo.dataset.variables['depth_0'].values.T
 DZ=nemo.dataset.variables['e3_0'].T
 DP=np.sum(DZ*Zd_mask,axis=2)
 nz=DZ.shape[2]
 if not 'density' in list(nemo.dataset.keys()):     
     nemo.construct_density(CT_AS=True,pot_dens=True)
 if not 'density_bar' in list(nemo.dataset.keys()):
     nemo.construct_density(CT_AS=True,rhobar=True,Zd_mask=Zd_mask.T,pot_dens=True)
 
 rho=nemo.dataset.variables['density'].values.T
 rhobar=nemo.dataset.variables['density_bar'].values.T
 rhobar_3d=np.repeat(rhobar[:,:,np.newaxis],nz,axis=2)
 PEA=np.sum(Z*(rho-rhobar_3d)*DZ*Zd_mask,axis=2)*g/DP
#%%
 return PEA
#%%
def calc_zdmask(nemo,nemo_w,Zmax):


 Z=nemo.dataset.variables['depth_0'].values.T
 ZW=nemo_w.dataset.variables['depth_0'].values.T
 mbot=nemo.dataset.variables['bottom_level'].values.T
 mask=mbot!=0
 ZZ=ZW[:,:,1:]
 ZZ[ZZ==0]=np.nan
 ZW[:,:,1:]=ZZ

 Zd_mask=np.zeros((Z.shape))
 nx,ny,nz=np.shape(Z)
 kmax=np.zeros((nx,ny)).astype(int)
 IIkmax=np.zeros(np.shape(Z))
#
#careful assumes mbot is 1st sea point above bed ie new definition
 for i in range(nx):
     for j in range(ny):
       if mask[i,j]==1:          
         Zd_mask[i,j,0:mbot[i,j]]=1 # mbot is not python style index so no +1
         kmax[i,j]=mbot[i,j]
         if ZW[i,j,mbot[i,j]]>Zmax:
           kkmax=np.max(np.where(ZW[i,j,:]<Zmax))
           Zd_mask[i,j,kkmax+1:]=0
           Zd_mask[i,j,kkmax]=(Zmax-ZW[i,j,kkmax])/(ZW[i,j,kkmax+1]-ZW[i,j,kkmax])
           kmax[i,j]=kkmax
           IIkmax[i,j,kkmax]=1
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
 
 
 
 

