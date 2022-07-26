#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 14:32:51 2022

@author: jholt
"""
import socket
isliv = 'livljobs' in socket.gethostname()
import numpy as np
import matplotlib.pylab as plt
import sys
import pandas as pd
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/work/Git/COAsT/')
 
import coast
import scipy.io
def error_stats(var,var_obs):
 I=np.nonzero((np.isfinite(var))&(np.isfinite(var_obs)))

 rmse=np.sqrt(np.mean((var[I]-var_obs[I])**2))
 me=np.mean(var[I]-var_obs[I])
 cost=rmse/np.std(var_obs[I])
 cc=np.corrcoef(var[I],var_obs[I])
 cc=cc[0,1]
 n=np.size(I)
 return rmse,me,cost,cc,n

J_offset=186 #account for extra rows in eORCA if data is made for normal ORCA
plot=False
metrics={}
A=np.load('/work/jholt/Git/DEV_jholt/GCO_Hydro/Data/LME_gridinfo_V4.npz')
a=scipy.io.loadmat('/projectsa/FASTNEt/jholt/HCIdata/ORCA025_ROAM_GLB_LMEmaskV4.mat')
LME_mask=a['LME_mask'][:,:].T
datapath_LME ='/work/jholt/Git/DEV_jholt/GCO_Hydro/Data/'
fn_bathymetry='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/eORCA025_bathy_meter.nc'
DATANAME='ORCA025'

nlme=66#len(A['i_min'])
Depth_lim=1000.
Assessdir='/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT/ORCA025-SE-NEMO/'
Outdir=Assessdir+'Errorstats/'
monrange=range(5,9)
lmonrange=len(monrange)
mons='NHS'

RUNNAMS=[
    'ORCA025-SE-NEMO_1980_1985_EXP_MES_NOTIDE',
    'ORCA025-SE-NEMO_1980_1985_EXP_MES_TIDE',
    'ORCA025-SE-NEMO_1980_1985_NOTIDE',
    'ORCA025-SE-NEMO_1980_1985_TIDE',
    'ORCA025-SE-NEMO_1979_1981_EXP_MES',
    'ORCA025-SE-NEMO_1979_1981_EXP_MES_WAV_NTM',
    'ORCA025-SE-NEMO_1979_1981_EXP_MES_WAV'
    ]


for RUNNAM in RUNNAMS:
    
    
    rmse_T=np.zeros(nlme)
    me_T=np.zeros(nlme)
    cost_T=np.zeros(nlme)
    cc_T=np.zeros(nlme)
    n_T=np.zeros(nlme)
    
    rmse_S=np.zeros(nlme)
    me_S=np.zeros(nlme)
    cost_S=np.zeros(nlme)
    cc_S=np.zeros(nlme)
    n_S=np.zeros(nlme)
    
    rmse_PEA=np.zeros(nlme)
    me_PEA=np.zeros(nlme)
    cost_PEA=np.zeros(nlme)
    cc_PEA=np.zeros(nlme)
    n_PEA=np.zeros(nlme)
    
    
    fn_nemo_dat=Assessdir +RUNNAM+'_SST_SSS_PEA_MonClimate.nc'
    f = coast.Gridded(fn_data= fn_nemo_dat)
    f_bathy=coast.Gridded(fn_data= fn_bathymetry)
    LMEs=[]
    for ilme in range(nlme):    
     LMENAM=A['DOMNAM'][ilme]
     i_min=A['i_min'][ilme]
     i_max=A['i_max'][ilme]
     j_min=A['j_min'][ilme]+J_offset
     j_max=A['j_max'][ilme]+J_offset
     j_min0=A['j_min'][ilme]
     j_max0=A['j_max'][ilme]
     Depth=f_bathy.dataset.variables['Bathymetry'].values[j_min:j_max+1,i_min:i_max+1]     
     PEAy=f.dataset.variables['PEAy'].values[monrange,j_min:j_max+1,i_min:i_max+1]
     SSTy=f.dataset.variables['SSTy'].values[monrange,j_min:j_max+1,i_min:i_max+1]
     SSSy=f.dataset.variables['SSSy'].values[monrange,j_min:j_max+1,i_min:i_max+1]
     mask=f.dataset.variables['bottom_level'].values[j_min:j_max+1,i_min:i_max+1] !=0
     mm=mask*LME_mask[j_min0:j_max0+1,i_min:i_max+1]==ilme+1
     mm=np.repeat(mm[np.newaxis,:,:],lmonrange,axis=0)
     Dmask=Depth<=Depth_lim
     Dmask=np.repeat(Dmask[np.newaxis,:,:],lmonrange,axis=0)
     mm=mm*Dmask
     SSSy[mm==0]=np.nan
     SSTy[mm==0]=np.nan
     PEAy[mm==0]=np.nan
    
     try:
    #%%
         EN4=np.load(datapath_LME + '/' + DATANAME  +'/'+ LMENAM + '_EN4mnthgrid_V3.npz') 
         SST_EN4=EN4['SST_g'].T[monrange,:,:]
         SSS_EN4=EN4['SSS_g'].T[monrange,:,:]
         PEA_EN4=EN4['PEA_g'].T[monrange,:,:]     
         SST_EN4[mm==0]=np.nan
         PEA_EN4[mm==0]=np.nan
         SSS_EN4[mm==0]=np.nan
         
         PEA_EN4[PEA_EN4<0]=0
         PEA_EN4[PEA_EN4>1e5]=np.nan     
    #     rmse_T[ilme],me_T[ilme],cost_T[ilme],cc_T[ilme],n_T[ilme]          =error_stats(SSTy.ravel(),SST_EN4.ravel())
    #     rmse_S[ilme],me_S[ilme],cost_S[ilme],cc_S[ilme],n_S[ilme]          =error_stats(SSSy.ravel(),SSS_EN4.ravel())
    #     rmse_PEA[ilme],me_PEA[ilme],cost_PEA[ilme],cc_PEA[ilme],n_PEA[ilme]=error_stats(PEAy.ravel(),PEA_EN4.ravel())
         metrics['rmse','SST',ilme],metrics['me','SST',ilme],metrics['cost','SST',ilme],metrics['cc','SST',ilme],metrics['N','SST',ilme]=error_stats(SSTy.ravel(),SST_EN4.ravel())
         metrics['rmse','SSS',ilme],metrics['me','SSS',ilme],metrics['cost','SSS',ilme],metrics['cc','SSS',ilme],metrics['N','SSS',ilme]=error_stats(SSSy.ravel(),SSS_EN4.ravel())
         metrics['rmse','PEA',ilme],metrics['me','PEA',ilme],metrics['cost','PEA',ilme],metrics['cc','PEA',ilme],metrics['N','PEA',ilme]=error_stats(PEAy.ravel(),PEA_EN4.ravel())
         LMEs=np.append(LMEs,ilme).astype(int)     
         if plot:     
             plt.figure(ilme)
             v1=PEAy[8,:,:]
             v2=PEA_EN4[8,:,:]
             vmax=np.nanmax(np.concatenate((v1.ravel(),v2.ravel())))
            
             plt.subplot(2,2,1)
             plt.pcolormesh(v1,vmax=vmax,vmin=0)
             plt.subplot(2,2,2)
             plt.pcolormesh(v2,vmax=vmax,vmin=0)
             plt.colorbar(orientation='vertical')
             
             v1=SSSy[8,:,:]
             v2=SSS_EN4[8,:,:]
             vmax=np.nanmax(np.concatenate((v1.ravel(),v2.ravel())))
             vmin=np.nanmin(np.concatenate((v1.ravel(),v2.ravel())))    
             plt.subplot(2,2,3)
             plt.pcolormesh(v1,vmax=vmax,vmin=vmin)
             plt.subplot(2,2,4)
             plt.pcolormesh(v2,vmax=vmax,vmin=vmin) 
             plt.colorbar(orientation='vertical')
             
    
    
         
    #%%     
     except:     
      print('failed:',ilme)
      
      
    #%%
      
      
    DD={}
    SS={}
    SN={}
    DD['LME']=[]
    vars=['SST','SSS','PEA']
    Metrics=['me','cost','cc','N']
    for var in vars:
       for Metric in Metrics:
        DD[var+' '+Metric]=[]
        SS[var+' '+Metric]=0
        SN[var+' '+Metric]=0
        
    for ilme in LMEs:
        DD['LME']=np.append(DD['LME'],A['DOMNAM'][ilme])
        for var in vars:
            for Metric in Metrics:
                DD[var+' '+Metric]=np.append(DD[var+' '+Metric],metrics[Metric,var,ilme])
                if np.isfinite(metrics[Metric,var,ilme]) and np.isfinite(metrics['N',var,ilme]):
                 SS[var+' '+Metric]=SS[var+' '+Metric]+metrics[Metric,var,ilme]*metrics['N',var,ilme]
                 SN[var+' '+Metric]=SN[var+' '+Metric]+metrics['N',var,ilme]
                 

    for var in vars:
        for Metric in Metrics:
            SS[var+' '+Metric]=SS[var+' '+Metric]/SN[var+' '+Metric]

            DD[var+' '+Metric]=np.append(DD[var+' '+Metric],SS[var+' '+Metric])
    DD['LME']=np.append(DD['LME'],'N Weighted Mean')                                     
    df=pd.DataFrame(DD)
    df=df.set_index('LME')
    df.to_csv(Outdir +'Errorstats_'+ RUNNAM +'_'+str(int(Depth_lim))+'m_'+mons+'_V1.csv')
    #writer = pd.ExcelWriter(Outdir +'Errorstats_'+ RUNNAM +'_V1.xlsx', engine='xlsxwriter')
    #df.to_excel(writer, sheet_name='Sheet1')
    #workbook  = writer.book
    #worksheet = writer.sheets['Sheet1']
    #writer.save()
                    
        
        
        
      