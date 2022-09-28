#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
THIS IS REDUNDENT USE make_SENEMO_climatology.py instead
Created on Thu Sep  5 16:27:35 2019

@author: jholt
"""
import socket
isliv = 'livljobs' in socket.gethostname()
import numpy as np
#This is not be needed if coast is installed propoerly
import sys
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/work/Git/COAsT/')
############    
# Needs SENEMO branch of coast
import coast
import xarray as xr

#Depth for PEA integration
Zmax=200
###############################################################################
def calc_pea(nemo,Zd_mask): #Redudent - can delete
# Calculates Potential Engergy Anomaly
#%%
 g=9.81
 Z=nemo.dataset.variables['depth_0'].values
 DZ=nemo.dataset.variables['e3_0'].values*Zd_mask
 DP=np.sum(DZ,axis=0) #water depth or Zmax , 

 nt=nemo.dataset.dims['t_dim']

 if not 'density' in list(nemo.dataset.keys()):     
        nemo.construct_density(CT_AS=True,pot_dens=True)
 if not 'density_bar' in list(nemo.dataset.keys()):
        nemo.construct_density(CT_AS=True,rhobar=True,Zd_mask=Zd_mask,pot_dens=True)
 rho=nemo.dataset.variables['density'].values #density 
 rho[np.isnan(rho)]=0
 rhobar=nemo.dataset.variables['density_bar'].values #density with depth-mean T and S

 z_axis=0
 if len(nemo.dataset['density'].shape) == 4:   # includes time as first axis
  Z=np.repeat(Z[np.newaxis,:,:,:],nt,axis=0)
  DZ=np.repeat(DZ[np.newaxis,:,:,:],nt,axis=0)
  DP=np.repeat(DP[np.newaxis,:,:],nt,axis=0)  

  z_axis=1
#%% 
 PEA=np.sum(Z*(rho-rhobar)*DZ,axis=z_axis)*g/DP

 return PEA
###############################################################################
def calc_zdmask(nemo_t,Zmax): #Redudent - can delete
#%%
 """
    Calculates a 3D mask so a specified level Zmax. 1 for sea; 0 for below sea bed
    and linearly ramped for last level
 """
 Z=nemo_t.dataset.variables['depth_0'].values
 e3t_0=nemo_t.dataset.variables['e3_0'].values
#calculate W-level - might want this done as stanbdard in gridded
 ZW = np.zeros_like(e3t_0)
 ZW[0, :, :] = 0.0
 ZW[1:, :, :] = np.cumsum(e3t_0, axis=0)[:-1, :, :]
 mbot=nemo_t.dataset.variables['bottom_level'].values.astype(int)
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
#%%
 return Zd_mask,kmax,Ikmax
###############################################################################         
class annual_climatology: #Redunded can delete
    def __init__ (self,gridded_t,gridded_t_out):
                
#%%
     Zd_mask,kmax,Ikmax=gridded_t.calculate_vertical_mask(Zmax) 
     #Calculate and save first time, otherwise read
     #try:     
     # A=np.load(domain_outpath + '/' +DOMNAM + '_' + EXPNAM + '_Zd_mask.npz')  
     # Zd_mask=A['Zd_mask']
     # kmax=A['kmax']
     # Ikmax=A['Ikmax']
     # print('read mask')
     #except:
     #   print('calculating mask')
     #   Zd_mask,kmax,Ikmax=gridded_t.calculate_vertical_mask(Zmax)
     #   np.savez(domain_outpath + '/' +DOMNAM + '_' + EXPNAM + '_Zd_mask.npz',Zd_mask=Zd_mask,kmax=kmax,Ikmax=Ikmax)
    #%%    
     ny=gridded_t.dataset.dims['y_dim']
     nx=gridded_t.dataset.dims['x_dim']
    
     nt=gridded_t.dataset.dims['t_dim']
    #%%
     print(nx,ny)
     SSTy=np.zeros((12,ny,nx))
     SSSy=np.zeros((12,ny,nx))
     PEAy=np.zeros((12,ny,nx))
     #NBTy=np.zeros((12,ny,nx))
     
    #%%  
     PEAy=np.zeros((12,ny,nx))
     nyear=int(nt/12)
     for iy in range(nyear):
      print('Calc PEA',iy)   
      it=np.arange((iy)*12,(iy)*12+12).astype(int)
      for im in range(12):
       itt=[it[im]]
       print(itt)
       gridded_t2=gridded_t.subset_as_copy(t_dim=itt)
       print('copied',im)
       PEA=coast.InternalTide(gridded_t2,gridded_t2)
       PEA.calc_pea(gridded_t2,Zd_mask)
       PEAy[im,:,:]=PEAy[im,:,:]+PEA.dataset['PEA'].values
     PEAy=PEAy/nyear 
    
    
     #need to find efficient method for bottom temperature
     #NBT=np.zeros((nt,ny,nx))
     #for it in range(nt): 
     #    NBT[it,:,:]=np.reshape(tmp[it,:,:,:].values.ravel()[Ikmax],(ny,nx))
     SST=gridded_t.dataset.variables['temperature'][:,0,:,:]
     SSS=gridded_t.dataset.variables['salinity'][:,0,:,:]    
    
     for im in range(12):
       print('Month',im)
       it=np.arange(im,nt,12).astype(int)
       SSTy[im,:,:]=np.mean(SST[it,:,:],axis=0)
       SSSy[im,:,:]=np.mean(SSS[it,:,:],axis=0)
     # NBTy[im,:,:]=np.mean(NBT[it,:,:],axis=0)
        # save hard work in netcdf file
     coords = {"Months":(("mon_dim"),np.arange(12).astype(int)),
                "latitude": (("y_dim", "x_dim"), gridded_t.dataset.latitude.values),
                "longitude": (("y_dim", "x_dim"), gridded_t.dataset.longitude.values),
            }
     dims = ["mon_dim","y_dim", "x_dim"]                
     attributes_SST = {"units": "o^C", "standard name": "Conservative Sea Surface Temperature"}
     attributes_SSS = {"units": "", "standard name": "Absolution Sea Surface Salinity"}
     attributes_PEA = {"units": "Jm^-3", "standard name": "Potential Energy Anomaly to 200m"}
     gridded_t_out.dataset['SSTy'] = xr.DataArray(np.squeeze(SSTy), coords=coords, dims=dims, attrs=attributes_SST) 
     gridded_t_out.dataset['SSSy'] = xr.DataArray(np.squeeze(SSSy), coords=coords, dims=dims, attrs=attributes_SSS) 
     gridded_t_out.dataset['PEAy'] = xr.DataArray(np.squeeze(PEAy), coords=coords, dims=dims, attrs=attributes_PEA) 
#%% 
##############################################################################
def NEMO_FileNames(dpath,runtype,ystart,ystop):  # redundent can delete
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
###############################################################################
if __name__ == '__main__':
    EXPNAMS=['EXP_SZT39_TAPER_TIDE','EXP_ZPS','EXP_SZT39_TAPER',
             'EXP_SZT39_TAPER_TKE','EXP_SZT51_NOTAPER','EXP_SZT39_NOTAPER']
    DOMCFGNAMS=['domain_cfg_ztaper_match.nc','domain_cfg_zps.nc',
                'domain_cfg_ztaper_match.nc','domain_cfg_ztaper_match.nc'
                ,'domain_cfg_51_noztaper_match_rmax15.nc','domain_cfg_noztaper_match.nc']
    
    EXPNAMS=['EXP_MES_TIDE','EXP_MES_NOTIDE','TIDE','NOTIDE']
    EXPNAMS=['EXP_MES',  'EXP_MES_WAV',  'EXP_MES_WAV_NTM'] 
    EXPNAMS=['EXP_MES_WAV_NTM']  
    #EXPNAMS=['EXP_SZT39_TAPER']
    for ik,EXPNAM in enumerate(EXPNAMS):
        print(EXPNAM)
        if isliv: #NOC-liverpool
            domain_datapath='/work/jholt/JASMIN//SENEMO/JDHA/' + EXPNAM +'/'
            domain_outpath='/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT/'
            
            if EXPNAM=='TIDE':
               domain_datapath='/work/jholt/JASMIN//SENEMO/TIDE/outputs/'

            if EXPNAM=='NOTIDE'    :
               domain_datapath='/work/jholt/JASMIN//SENEMO/NOTIDE/'
                           
        else: #JASMIN
            domain_datapath='/gws/nopw/j04/class_vol2/senemo/jdha/SHORT_TESTS/' + EXPNAM +'/'
            fn_nemo_dom='/gws/nopw/j04/class_vol2/senemo/jdha/SHORT_TESTS/EXP_MES/domain_cfg_r015-r010_007_004v2.nc'
            if EXPNAM=='TIDE':
               domain_datapath='/gws/nopw/j04/class_vol2/senemo/dbyrne/EXP_REF_TIDE/outputs/'
               fn_nemo_dom='/gws/nopw/j04/class_vol2/senemo/cwilso01/senemo/EXP_REF_NOTIDE/domcfg_eORCA025_v2.nc'
 
            if EXPNAM=='NOTIDE'    :
               domain_datapath='/gws/nopw/j04/class_vol2/senemo/cwilso01/senemo/EXP_REF_NOTIDE/means/monthly/'
               fn_nemo_dom='/gws/nopw/j04/class_vol2/senemo/cwilso01/senemo/EXP_REF_NOTIDE/domcfg_eORCA025_v2.nc'

            domain_outpath='/home/users/jholt/work/SENEMO/ASSESSMENT/'
        domain_path=domain_datapath
        
        
        
        #fn_nemo_dom=domain_path+DOMCFGNAMS[ik]
        #fn_nemo_dom='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc'
        
        print(fn_nemo_dom)
        
        #make list of filenames
        ystart=1979
        ystop=1979
        fn_nemo_dat= coast.nemo_filenames(domain_datapath,'SENEMO',ystart,ystop)            

        fn_config_t_grid='../Config/senemo_grid_t.json'    
        
        #input datasets
        nemo = coast.Gridded(fn_data= fn_nemo_dat, fn_domain = fn_nemo_dom, config=fn_config_t_grid,multiple=True)
<<<<<<< HEAD
        nemo = nemo. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))  
=======
        #nemo = nemo. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))
        #fix for nasty bug 
        nemo_dom=coast.Gridded(fn_domain = fn_nemo_dom, config=fn_config_t_grid)
        #nemo_dom = nemo_dom. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))        
        nemo.dataset['e3_0']=nemo_dom.dataset['e3_0']
>>>>>>> fbd186d847080e062b2c7d5edc2b69c4cc56554b
        #Place to output data
        nemo_out=coast.Gridded(fn_domain = fn_nemo_dom, config=fn_config_t_grid)
        #nemo_out = nemo_out. subset_as_copy(y_dim=range(86,1000),x_dim=range(1080,1180))  
        fn_nameout=EXPNAM+ '_SST_SSS_PEA_MonClimate.nc'
        DOMNAM='ORCA025-SE-NEMO'
        
        fn_out=domain_outpath+'/'+DOMNAM +'/' +DOMNAM+'_1979_1979_'+fn_nameout
        print('running')
        #%% do the hard work
        coast.Annual_Climatology(nemo,nemo_out,Zmax=200)        
        nemo_out.dataset.to_netcdf(fn_out)
    
    
    
    
