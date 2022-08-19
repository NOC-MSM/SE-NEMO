#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 16:46:08 2022

@author: jholt
"""
import socket
isliv = 'livljobs' in socket.gethostname()
if isliv:
    coast_dir='/login/jholt/work/Git/COAsT/'
else:
    coast_dir='/home/users/jholt/work/Git/COAsT/'
 
import sys
sys.path.insert(0,coast_dir)
import matplotlib.pylab as plt
import coast
import numpy as np
import xarray as xr
import scipy.io
fn_config_t_grid='./example_nemo_grid_t.json'
fn_config_f_grid='./example_nemo_grid_f.json'
fn_config_u_grid='./example_nemo_grid_u.json'
fn_config_v_grid='./example_nemo_grid_v.json'


names,dpaths,DOMS,_  = coast. experiments(experiments='experiments1.json')  

def flx_contour(nemo_f,nemo_u,nemo_v):
    contours, no_contours = coast.Contour.get_contours(nemo_f, 300)
    y_ind, x_ind, contour = coast.Contour.get_contour_segment(nemo_f, contours[0], [57.4, -9.6], [61.64, 1.79])
    cont_f = coast.ContourF(nemo_f, y_ind, x_ind, 300)
    cont_f.calc_cross_contour_flow(nemo_u, nemo_v)

    Un=cont_f.data_cross_flow.normal_velocities.values
 #%%   
    Un[Un==0]=np.nan
    if len(Un.shape)==3:
        Untm=np.nanmean(Un,axis=0)
        Unmean=np.nanmean(Untm,axis=1)
    else:
        Unmean=np.nanmean(Un,axis=1)            
    Z=-cont_f.data_cross_flow.depth_0.values
    Zmean=np.nanmean(Z,axis=1)
    r_dim=cont_f.data_cross_flow.r_dim.values
    return Un,Unmean,Z,Zmean
#%%
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
def lightcolormap(Np,nt):
    
 cmap0=cm.get_cmap('BrBG_r',lut=Np+nt*2)
 colors = cmap0(np.arange(cmap0.N))
 colors1=colors[nt:cmap0.N-nt]
 cmap1=LinearSegmentedColormap.from_list('cmap1',colors1,cmap0.N-nt*2)
 return cmap1
#%%
def plot_surface_circulation(nemo_u,nemo_v,nemo_t,mask,name, co_located=False,Vmax=0.16,Np=3): 
#%%
 nx=nemo_u.dataset.x_dim.size
 ny=nemo_u.dataset.y_dim.size
 VS=np.zeros((ny,nx))
 US=np.zeros((ny,nx))


 if len(nemo_v.dataset.v_velocity.dims)==2:
     vs=nemo_v.dataset.v_velocity.values
     us=nemo_u.dataset.u_velocity.values
 elif len(nemo_v.dataset.v_velocity.dims)==3:
     vs=nemo_v.dataset.v_velocity[:,:,:].mean(dim="t_dim").values
     us=nemo_u.dataset.u_velocity[:,:,:].mean(dim="t_dim").values
 else:     
     vs=nemo_v.dataset.v_velocity[:,0,:,:].mean(dim="t_dim")
     us=nemo_u.dataset.u_velocity[:,0,:,:].mean(dim="t_dim")
 
 if co_located:
     VS=vs
     US=us
 else:
     VS[:-1,:]=0.5*(vs[:-1,:]+vs[1:,:])
     US[:,:-1]=0.5*(us[:,:-1]+us[:,1:]) 

 SP=np.sqrt(US**2+VS**2)
 US[SP<0.02]=np.nan
 VS[SP<0.02]=np.nan
 US=US/SP
 VS=VS/SP
 


 p=np.ma.masked_where(mask==0,SP)
 u=np.ma.masked_where(mask[0::Np,0::Np]==0,US[0::Np,0::Np])
 v=np.ma.masked_where(mask[0::Np,0::Np]==0,VS[0::Np,0::Np])
 #x=nemo_t.dataset.longitude[0::Np,0::Np]
 #y=nemo_t.dataset.latitude[0::Np,0::Np]
 X=nemo_t.dataset.longitude
 Y=nemo_t.dataset.latitude

 cmap1=lightcolormap(int(Vmax*100),2)
 cmap1.set_bad([0.75,0.75,0.75])
 fig=plt.figure(); fig.clear()
 plt.pcolormesh(p,cmap=cmap1)
 #plt.pcolormesh(X,Y,p,cmap=cmap1)
 x=np.arange(0,nx,Np)
 y=np.arange(0,ny,Np)
 
 plt.clim([0,Vmax])
 plt.colorbar(orientation='vertical',cmap=cmap1)
 #plt.quiver(x,y,u,v,color=[0.1,0.1,0.1],headwidth=4,scale=50)
 plt.quiver(x,y,u,v,color=[0.1,0.1,0.1],headwidth=4,scale=50)
 plt.title('Surface Currents ' + name)




#%%
def plot_flx_contour(Unmean,Zmean,names):
    plt.plot(Unmean[names[0]],Zmean[names[0]],'o-',
             Unmean[names[1]],Zmean[names[1]],'o-',
             Unmean[names[2]],Zmean[names[2]],'o-',
             Unmean[names[3]],Zmean[names[3]],'o-')       
             
    plt.ylim((-350,0))
    plt.legend(names)
######################################################################
if __name__ == '__main__':
    Unmean={}
    Untm={}
    Zmean={}
    Un={}
    Z={}
    ystart=1993
    ystop=2019
    A=np.load('../Data/LME_gridinfo_V4.npz')
    a=scipy.io.loadmat('../Data/ORCA025_ROAM_GLB_LMEmaskV4.mat')
    nlme=66
    J_offset=186 #account for extra rows in eORCA if data is made for normal ORCA
    LME_mask=a['LME_mask'][:,:].T
    lmelist=np.array([14])-1

    SEASON='JAS'
    YEARS="{0}_{1}".format(ystart,ystop)
    plt.close('all')
    ny=ystop-ystart+1
    isea=[6,7,8]
    iseason=np.array([])
    for iss in isea:        
        iseason=np.append(iseason,np.arange(iss,12*ny,12))
    iseason=(np.sort(iseason)).astype(int)                        
    for i,name in enumerate(names): 
        print(name)
        print(dpaths[i])
        fn_nemo_dat_u= coast.nemo_filenames(dpaths[i],'SENEMO',ystart,ystop,grid='U') 
        fn_nemo_dat_v= coast.nemo_filenames(dpaths[i],'SENEMO',ystart,ystop,grid='V') 
        fn_nemo_dat_u=    fn_nemo_dat_u[iseason]       
        fn_nemo_dat_v=    fn_nemo_dat_v[iseason]   
        fn_nemo_dom=DOMS[i]
    #    fn_nemo_dat_u=dpaths[i]+'SENEMO_1m_19800101_19801231_grid_U_198007-198007.nc'
    #    fn_nemo_dat_v=dpaths[i]+'SENEMO_1m_19800101_19801231_grid_V_198007-198007.nc'
    
        
        
        
        nemo_t = coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_t_grid)#,calc_bathy=True)    
        nemo_f = coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_f_grid)#,calc_bathy=True)
        nemo_u = coast.Gridded(fn_data=fn_nemo_dat_u, fn_domain=fn_nemo_dom, config=fn_config_u_grid,multiple=True)
        nemo_v = coast.Gridded(fn_data=fn_nemo_dat_v, fn_domain=fn_nemo_dom, config=fn_config_v_grid,multiple=True)
        for ilme in lmelist:    
            LMENAM=A['DOMNAM'][ilme]
            imin=A['i_min'][ilme]
            imax=A['i_max'][ilme]
            jmin=A['j_min'][ilme]+J_offset
            jmax=A['j_max'][ilme]+J_offset
            jmin0=A['j_min'][ilme]
            jmax0=A['j_max'][ilme]
            print(LMENAM)
            
            jmin,jmax=[860,1015]
            imin,imax=[1080,1180]
            REGION=LMENAM
            REGION='NWS'            
            nemo_f1=nemo_f.subset_as_copy(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
            nemo_u1=nemo_u.subset_as_copy(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
            nemo_v1=nemo_v.subset_as_copy(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
            nemo_t1=nemo_t.subset_as_copy(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
            mask=nemo_t1.dataset.bathymetry != 0
            Name=name+' '+SEASON+' '+YEARS+' '+REGION 
            plot_surface_circulation(nemo_u1,nemo_v1,nemo_t1,mask,Name)
            plt.savefig('../Figures/Circulation/Surface_Currents_' + Name.replace(' ','_')+'.png')
    #    Un[name],Unmean[name],Z[name],Zmean[name]=flx_contour(nemo_f,nemo_u,nemo_v)
    
    
    
    
    
    #plot_flx_contour(Unmean,Zmean,names)
        
