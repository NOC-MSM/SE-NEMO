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
fn_config_t_grid=coast_dir+'/config/example_nemo_grid_t.json'
fn_config_f_grid=coast_dir+'/config/example_nemo_grid_f.json'
fn_config_u_grid=coast_dir+'/config/example_nemo_grid_u.json'
fn_config_v_grid=coast_dir+'/config/example_nemo_grid_v.json'

names,dpaths,DOMS,_  = coast. experiments(experiments='experiments.json')  

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
def plot_surface_circulation(nemo_u,nemo_v,nemo_t,name): 

 nx=nemo_u.dataset.x_dim.size
 ny=nemo_u.dataset.y_dim.size
 VS=np.zeros((ny,nx))
 US=np.zeros((ny,nx))


 vs=nemo_v.dataset.v_velocity[:,0,:,:].mean(dim="t_dim")
 us=nemo_u.dataset.u_velocity[:,0,:,:].mean(dim="t_dim")
 VS[:-1,:]=0.5*(vs[:-1,:]+vs[1:,:])
 US[:,:-1]=0.5*(us[:,:-1]+us[:,1:]) 
 SP=np.sqrt(US**2+VS**2)
 US[SP<0.02]=np.nan
 VS[SP<0.02]=np.nan
 US=US/SP
 VS=VS/SP
 Np=3

 D=nemo_t.dataset.bathymetry
 p=np.ma.masked_where(D==0,SP)
 u=np.ma.masked_where(D[0::Np,0::Np]==0,US[0::Np,0::Np])
 v=np.ma.masked_where(D[0::Np,0::Np]==0,VS[0::Np,0::Np])
 x=nemo_t.dataset.longitude[0::Np,0::Np]
 y=nemo_t.dataset.latitude[0::Np,0::Np]
 X=nemo_t.dataset.longitude
 Y=nemo_t.dataset.latitude

 cmap1=lightcolormap(16,2)
 cmap1.set_bad([0.75,0.75,0.75])
 fig=plt.figure(); fig.clear()
 plt.pcolormesh(p,cmap=cmap1)
 #plt.pcolormesh(X,Y,p,cmap=cmap1)
 I=np.arange(0,nx,Np)
 J=np.arange(0,ny,Np)
 
 plt.clim([0,0.16])
 plt.colorbar(orientation='vertical',cmap=cmap1)
 #plt.quiver(x,y,u,v,color=[0.1,0.1,0.1],headwidth=4,scale=50)
 plt.quiver(I,J,u,v,color=[0.1,0.1,0.1],headwidth=4,scale=50)
 plt.title('Surface Currents ' + name)
 plt.savefig('../Figures/Circulation/Surface_Currents_' + name.replace(' ','_')+'.png')



#%%
def plot_flx_contour(Unmean,Zmean,names):
    plt.plot(Unmean[names[0]],Zmean[names[0]],'o-',
             Unmean[names[1]],Zmean[names[1]],'o-',
             Unmean[names[2]],Zmean[names[2]],'o-',
             Unmean[names[3]],Zmean[names[3]],'o-')       
             
    plt.ylim((-350,0))
    plt.legend(names)
######################################################################
Unmean={}
Untm={}
Zmean={}
Un={}
Z={}
ystart=1981
ystop=1981
jmin,jmax=[860,1010]
imin,imax=[1080,1180]
plt.close('all')
for i,name in enumerate(names): 
    print(name)
    print(dpaths[i])
    fn_nemo_dat_u= coast.nemo_filenames(dpaths[i],'SENEMO',ystart,ystop,grid='U') 
    fn_nemo_dat_v= coast.nemo_filenames(dpaths[i],'SENEMO',ystart,ystop,grid='V') 
    fn_nemo_dat_u=    fn_nemo_dat_u[6:9]       
    fn_nemo_dat_v=    fn_nemo_dat_v[6:9]   
    fn_nemo_dom=DOMS[i]
#    fn_nemo_dat_u=dpaths[i]+'SENEMO_1m_19800101_19801231_grid_U_198007-198007.nc'
#    fn_nemo_dat_v=dpaths[i]+'SENEMO_1m_19800101_19801231_grid_V_198007-198007.nc'

    
    
    
    nemo_t = coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_t_grid)#,calc_bathy=True)    
    nemo_f = coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_f_grid)#,calc_bathy=True)
    nemo_u = coast.Gridded(fn_data=fn_nemo_dat_u, fn_domain=fn_nemo_dom, config=fn_config_u_grid,multiple=True)
    nemo_v = coast.Gridded(fn_data=fn_nemo_dat_v, fn_domain=fn_nemo_dom, config=fn_config_v_grid,multiple=True)
    
    nemo_f.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    nemo_u.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    nemo_v.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    nemo_t.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
    plot_surface_circulation(nemo_u,nemo_v,nemo_t,name+' JAS 1981')
#    Un[name],Unmean[name],Z[name],Zmean[name]=flx_contour(nemo_f,nemo_u,nemo_v)





#plot_flx_contour(Unmean,Zmean,names)
    
