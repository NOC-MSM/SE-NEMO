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
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
 
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

A=np.load('../Data/LME_gridinfo_V4.npz')
a=scipy.io.loadmat('../Data/ORCA025_ROAM_GLB_LMEmaskV4.mat')
LME_mask=a['LME_mask'][:,:].T
if isliv:
 datapath_LME ='/work/jholt/Git/DEV_jholt/GCO_Hydro/Data/'
 fn_bathymetry='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/eORCA025_bathy_meter.nc'
else:
 datapath_LME ='/home/users/jholt/work/SENEMO/ASSESSMENT/EN4/'
 fn_bathymetry='/home/users/jholt/work/SENEMO/senemo/INPUTS/eORCA025_bathy_meter.nc'    
DATANAME='ORCA025'

nlmes=66#len(A['i_min'])
Depth_lim=400.
if isliv:
 Assessdir='/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT/ORCA025-SE-NEMO/' 
else:
 Assessdir='/home/users/jholt//work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/'
Outdir=Assessdir+'Errorstats/'
monrange=range(0,12)

lmonrange=len(monrange)
mons='NHS'
mons='YR'

RUNNAMS=[
'ORCA025-SE-NEMO_1990_2019_EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2',
'ORCA025-SE-NEMO_1990_2019_ZPS_TIDE',
'ORCA025-SE-NEMO_1990_2019_ZPS_NOTIDE'  
    ]

SS={}
SN={}
LME_Clusters='../Data/LME_Clusters_eORCA025.csv'
clusters = pd.read_csv(LME_Clusters)
nclusters=clusters.values.shape[0]
metrics={}
for Vars in ['SST','SSS','PEA']:
    for Mets in ['rmse','me','cost','cc','N']:
      metrics[Mets,Vars]=np.zeros((nclusters,len(RUNNAMS))) *np.nan          
for iR,RUNNAM in enumerate(RUNNAMS):
#%%    
    
    
    fn_nemo_dat=Assessdir +RUNNAM+'_SST_SSS_PEA_MonClimate.nc'
    f = coast.Gridded(fn_data= fn_nemo_dat)
    f_bathy=coast.Gridded(fn_data= fn_bathymetry)
    LMEs=[]
    vnames=['PEAy','SSTy','SSSy']
    #if iR==11:
    vnames=['PEA_monthy_clim', 'SST_monthy_clim','SSS_monthy_clim']
    iwant_clusters=np.concatenate((np.arange(1,17),np.arange(18,23)))
    for icluster in iwant_clusters:
#%%
      #try:  
         lmes=clusters.values[icluster,10:]
         ilmes=np.array([])
         for LME in lmes:
             try:
                 ilmes=np.append(ilmes,int(LME[-2:]))
             except:
                 ilmes=np.append(ilmes,-1)
                 
     #%% 
           
         Cluster_Name=clusters.values[icluster,1]
         #eORCA025 index of cluster
         i_min1=int(clusters.values[icluster,2])
         i_max1=int(clusters.values[icluster,3])
         j_min1=int(clusters.values[icluster,4])
         j_max1=int(clusters.values[icluster,5])
         
    
         Depth=f_bathy.dataset.variables['Bathymetry'].values[j_min1:j_max1+1,i_min1:i_max1+1]     
         PEAy=f.dataset.variables[vnames[0]].values[monrange,j_min1:j_max1+1,i_min1:i_max1+1]
         SSTy=f.dataset.variables[vnames[1]].values[monrange,j_min1:j_max1+1,i_min1:i_max1+1]
         SSSy=f.dataset.variables[vnames[2]].values[monrange,j_min1:j_max1+1,i_min1:i_max1+1]
         mask=f.dataset.variables['bottom_level'].values[j_min1:j_max1+1,i_min1:i_max1+1] !=0
         #mm=mask*LME_mask[j_min0:j_max0+1,i_min:i_max+1]==ilme+1
         #mm=np.repeat(mm[np.newaxis,:,:],lmonrange,axis=0)
         Dmask=Depth<=Depth_lim
         Dmask=np.repeat(Dmask[np.newaxis,:,:],lmonrange,axis=0)
         mm=Dmask
         SSSy[mm==0]=np.nan
         SSTy[mm==0]=np.nan
         PEAy[mm==0]=np.nan
         SST_EN4=np.zeros(PEAy.shape)
         SSS_EN4=np.zeros(PEAy.shape)
         PEA_EN4=np.zeros(PEAy.shape)     
         mask=np.zeros((PEAy.shape[1],PEAy.shape[2]))
         for ilme,LMENAM in enumerate(lmes):
            if isinstance(LMENAM,str):
                
                
                
                LMEmask=LME_mask[(j_min1-J_offset):(j_max1-J_offset+1),i_min1:i_max1+1]
                EN4=coast.Gridded(datapath_LME + '/' + DATANAME  +'/'+ LMENAM + '_EN4mnthgrid_V3.nc')
                #ORCA025 index of LME bounds
                lims=EN4.dataset.lims.values
                #conver to eORCA grid
                lims[[2,3]]=lims[[2,3]]+J_offset
                
                lims1=np.copy(lims)
                lims1[0]=np.max([lims[0],i_min1])
                lims1[1]=np.min([lims[1],i_max1])
                lims1[2]=np.max([lims[2],j_min1])
                lims1[3]=np.min([lims[3],j_max1])                
                #index of the LME data
                i_min2=lims1[0]-lims[0]
                i_max2=lims1[1]-lims[0]
                j_min2=lims1[2]-lims[2]
                j_max2=lims1[3]-lims[2]                
                #Index of LME portion of CLuster
                i_min=lims1[0]-i_min1
                i_max=lims1[1]-i_min1
                j_min=lims1[2]-j_min1
                j_max=lims1[3]-j_min1
                mask[LMEmask==ilmes[ilme]]=1
                SST_EN4_LME=EN4.dataset.SST_g.values.T[monrange,:,:]
                SSS_EN4_LME=EN4.dataset.SSS_g.values.T[monrange,:,:]
                PEA_EN4_LME=EN4.dataset.PEA_g.values.T[monrange,:,:]            
                SST_EN4_LME[~np.isfinite(SST_EN4_LME)]=0.
                SSS_EN4_LME[~np.isfinite(SSS_EN4_LME)]=0.
                PEA_EN4_LME[~np.isfinite(PEA_EN4_LME)]=0.
                
                #nx=SST_EN4_LME.shape[2]
                #i1=0
                #i2=nx
                #if lims[0]<i_min1:
                #    i1=i_min1-lims[0]
                #    i_min=0
                # 
                #if lims[1]>i_max1:
                #    i2=nx-(lims[1]-i_max1+1)
                
                mm=LMEmask[j_min:j_max+1,i_min:i_max+1]==ilmes[ilme]

                SST_EN4[:,j_min:j_max+1,i_min:i_max+1]=(SST_EN4[:,j_min:j_max+1,i_min:i_max+1]
                +SST_EN4_LME[:,j_min2:j_max2+1,i_min2:i_max2+1]*mm)
                SSS_EN4[:,j_min:j_max+1,i_min:i_max+1]=(SSS_EN4[:,j_min:j_max+1,i_min:i_max+1]
                +SSS_EN4_LME[:,j_min2:j_max2+1,i_min2:i_max2+1]*mm)
                PEA_EN4[:,j_min:j_max+1,i_min:i_max+1]=(PEA_EN4[:,j_min:j_max+1,i_min:i_max+1]
                +PEA_EN4_LME[:,j_min2:j_max2+1,i_min2:i_max2+1]*mm)
                #        EN4=np.load(datapath_LME + '/' + DATANAME  +'/'+ LMENAM + '_EN4mnthgrid_V3.npz')         
             #        SST_EN4=EN4['SST_g'].T[monrange,:,:]
             #        SSS_EN4=EN4['SSS_g'].T[monrange,:,:]
             #        PEA_EN4=EN4['PEA_g'].T[monrange,:,:]     
            #         SST_EN4[mm==0]=np.nan
            #         PEA_EN4[mm==0]=np.nan
            #         SSS_EN4[mm==0]=np.na
            #except:     
            #        print('failed:',LMENAM)
         
         
         mask=np.repeat(mask[np.newaxis,:,:],12,axis=0)
        
         SST_EN4[SST_EN4==0]=np.nan
         SSS_EN4[SSS_EN4==0]=np.nan
         PEA_EN4[PEA_EN4==0]=np.nan
         
         SST_EN4[mask==0]=np.nan
         SSS_EN4[mask==0]=np.nan
         PEA_EN4[mask==0]=np.nan
         
         
         PEA_EN4[PEA_EN4<0]=0
         PEA_EN4[PEA_EN4>1e5]=np.nan     
        #     rmse_T[ilme],me_T[ilme],cost_T[ilme],cc_T[ilme],n_T[ilme]          =error_stats(SSTy.ravel(),SST_EN4.ravel())
        #     rmse_S[ilme],me_S[ilme],cost_S[ilme],cc_S[ilme],n_S[ilme]          =error_stats(SSSy.ravel(),SSS_EN4.ravel())
        #     rmse_PEA[ilme],me_PEA[ilme],cost_PEA[ilme],cc_PEA[ilme],n_PEA[ilme]=error_stats(PEAy.ravel(),PEA_EN4.ravel())
         metrics['rmse','SST'][icluster,iR],metrics['me','SST'][icluster,iR],metrics['cost','SST'][icluster,iR],metrics['cc','SST'][icluster,iR],metrics['N','SST'][icluster,iR]=error_stats(SSTy.ravel(),SST_EN4.ravel())
         metrics['rmse','SSS'][icluster,iR],metrics['me','SSS'][icluster,iR],metrics['cost','SSS'][icluster,iR],metrics['cc','SSS'][icluster,iR],metrics['N','SSS'][icluster,iR]=error_stats(SSSy.ravel(),SSS_EN4.ravel())
         metrics['rmse','PEA'][icluster,iR],metrics['me','PEA'][icluster,iR],metrics['cost','PEA'][icluster,iR],metrics['cc','PEA'][icluster,iR],metrics['N','PEA'][icluster,iR]=error_stats(PEAy.ravel(),PEA_EN4.ravel())
         LMEs=np.append(LMEs,icluster).astype(int)     
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
             
      #except:
      #   print('failed',icluster,clusters.values[icluster,1])
    #%%
      
    #Table for each experiment  
    DD={}

    DD['LME']=[]
    vars=['SST','SSS','PEA']
    Metrics=['me','cost','cc','N']
    for var in vars:
       for Metric in Metrics:
        DD[var+' '+Metric]=[]
        SS[var+' '+Metric,iR]=0
        SN[var+' '+Metric,iR]=0
        
    for icluster in range(nclusters):
        DD['LME']=np.append(DD['LME'],A['DOMNAM'][icluster])
        for var in vars:
            for Metric in Metrics:
                DD[var+' '+Metric]=np.append(DD[var+' '+Metric],metrics[Metric,var][icluster,iR])
                if np.isfinite(metrics[Metric,var][icluster,iR]) and np.isfinite(metrics['N',var][icluster,iR]):
                 SS[var+' '+Metric,iR]=SS[var+' '+Metric,iR]+metrics[Metric,var][icluster,iR]*metrics['N',var][icluster,iR]
                 SN[var+' '+Metric,iR]=SN[var+' '+Metric,iR]+metrics['N',var][icluster,iR]
                 

    for var in vars:
        for Metric in Metrics:
            SS[var+' '+Metric,iR]=SS[var+' '+Metric,iR]/SN[var+' '+Metric,iR]

            DD[var+' '+Metric]=np.append(DD[var+' '+Metric],SS[var+' '+Metric,iR])
            
    DD['LME']=np.append(DD['LME'],'N Weighted Mean')                                     
    df=pd.DataFrame(DD)
    df=df.set_index('LME')
    df.to_csv(Outdir +'Errorstats_clusters_'+ RUNNAM +'_'+str(int(Depth_lim))+'m_'+mons+'_V2.csv')
    #writer = pd.ExcelWriter(Outdir +'Errorstats_'+ RUNNAM +'_V1.xlsx', engine='xlsxwriter')
    #df.to_excel(writer, sheet_name='Sheet1')
    #workbook  = writer.book
    #worksheet = writer.sheets['Sheet1']
    #writer.save()
#%%
plt.close('all')
runs=['GS1p2_full','GS1p1_tide','GS1p0']

from matplotlib import cm    
cmap0=cm.get_cmap('BrBG_r',lut=16)
#vmin=0
#vmax=1.5
#vmin=-50
#vmax=50
#Metric='me'
#Var='PEA'
vmin=[ 1.5,0  ,0.6, 1.5,0  ,0.6,-80,0  ,0.6]
vmax=[-1.5,1.5,1.0,-1.5,1.5,1.0, 80,1.5,1.0]
isp=0
for Var in vars:    
    for Metric in Metrics[:3]:
        plt.figure()
        plt.pcolormesh(metrics[Metric,Var],cmap=cmap0,vmin=vmin[isp],vmax=vmax[isp])
        plt.yticks(ticks=np.arange(nclusters)+0.5,labels=clusters.values[:,1])
        plt.xticks(ticks=np.arange(3)+0.5,labels=runs)
        
        plt.title(Var + ' ' +Metric)
        ax=plt.gca()   
        ax.set_position([0.25, 0.11, 0.6, 0.8])
        cbax=plt.axes()
        cbax.set_position([0.87, 0.11, 0.04, 0.8])
        plt.colorbar(cax=cbax,orientation='vertical')
        isp=isp+1
        fname='../Figures/{0}_{1}_Clim.png'.format(Var,Metric)
        plt.savefig(fname)
#%%                    
    #Table for each LME  
for ilme in LMEs:
    LMENAME=A['DOMNAM'][ilme]
    DD={}

    DD['RUNNAM']=[]
    vars=['SST','SSS','PEA']
    Metrics=['me','cost','cc','N']
    for var in vars:
       for Metric in Metrics:
        DD[var+' '+Metric]=[]

        
    for iR,RUNNAM in enumerate(RUNNAMS):
        
        DD['RUNNAM']=np.append(DD['RUNNAM'],RUNNAM)
        for var in vars:
            for Metric in Metrics:
                DD[var+' '+Metric]=np.append(DD[var+' '+Metric],metrics[Metric,var,ilme,iR])
             
    df=pd.DataFrame(DD)
    df=df.set_index('RUNNAM')
    df.to_csv(Outdir +'Errorstats_'+ LMENAME +'_'+str(int(Depth_lim))+'m_'+mons+'_V4.csv')        
#%% Summary
DD={}
DD['RUNNAM']=[]
vars=['SST','SSS','PEA']
Metrics=['me','cost','cc','N']
for var in vars:
   for Metric in Metrics:
    DD[var+' '+Metric]=[]
for iR,RUNNAM in enumerate(RUNNAMS):    
    DD['RUNNAM']=np.append(DD['RUNNAM'],RUNNAM)        
    for var in vars:
            for Metric in Metrics:
                DD[var+' '+Metric]=np.append(DD[var+' '+Metric],SS[var+' '+Metric,iR]) 
df=pd.DataFrame(DD)
df=df.set_index('RUNNAM')
df.to_csv(Outdir +'Errorstats_clusters_'+ 'N_Weight_Mean' +'_'+str(int(Depth_lim))+'m_'+mons+'_V4.csv')
#%%
VAR_ERR=np.zeros((len(RUNNAMS),nclusters))*np.nan
for iR,RUNNAM in enumerate(RUNNAMS):
    for ilme in range(nclusters):
        try:
            VAR_ERR[iR,ilme]=metrics['cost','PEA',ilme,iR]
        except:
            continue





        