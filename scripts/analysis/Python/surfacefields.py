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
def mean_surface_circulation(nemo_u,nemo_v,nemo_t,mask,
                             co_located=False): 
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
     VS[1:,:]=0.5*(vs[:-1,:]+vs[1:,:])
     US[:,1:]=0.5*(us[:,:-1]+us[:,1:]) 

 SP=np.sqrt(US**2+VS**2)
 US[SP<0.02]=np.nan
 VS[SP<0.02]=np.nan
 US=US/SP
 VS=VS/SP
 return SP, US, VS

def mean_surface_TS(nemo_t,mask): 
#%%
 nx=nemo_t.dataset.x_dim.size
 ny=nemo_t.dataset.y_dim.size


 if len(nemo_t.dataset.temperature.dims)==2:
     TMP=nemo_t.dataset.temperature.values
     SAL=nemo_t.dataset.salinity.values
 elif len(nemo_t.dataset.temperature.dims)==3:
     TMP=nemo_t.dataset.temperature[:,:,:].mean(dim="t_dim").values
     SAL=nemo_t.dataset.salinity[:,:,:].mean(dim="t_dim").values
 else:     
     TMP=nemo_t.dataset.temperature[:,0,:,:].mean(dim="t_dim").values
     SAL=nemo_t.dataset.salinity[:,0,:,:].mean(dim="t_dim").values

 return TMP,SAL

#%%
def plot_surface_circulation(SP, US, VS,nemo_t,mask,name
                             ,Vmax=0.16,Np=3
                             ,headwidth=4,scale=50                            
                             ,**kwargs): 
 nx=nemo_t.dataset.x_dim.size
 ny=nemo_t.dataset.y_dim.size   
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
 plt.quiver(x,y,u,v,color=[0.1,0.1,0.1],headwidth=headwidth,scale=scale,**kwargs)
 plt.title('Surface Currents ' + name)
#%% 
def plot_surface_T_field(Var,nemo_t,mask,name,Vmin,Vmax):
    nx=nemo_t.dataset.x_dim.size
    ny=nemo_t.dataset.y_dim.size   
    Var=np.ma.masked_where(mask==0,Var)
    cmap1=lightcolormap(int(Vmax*100),2)
    cmap1.set_bad([0.75,0.75,0.75])
    fig=plt.figure(); fig.clear()
    plt.pcolormesh(Var,cmap=cmap1,vmin=Vmin,vmax=Vmax)
    plt.colorbar(orientation='vertical',cmap=cmap1)
    plt.title('Surface Sal Anom ' + name)


#%%
def plot_flx_contour(Unmean,Zmean,names):
    plt.plot(Unmean[names[0]],Zmean[names[0]],'o-',
             Unmean[names[1]],Zmean[names[1]],'o-',
             Unmean[names[2]],Zmean[names[2]],'o-',
             Unmean[names[3]],Zmean[names[3]],'o-')       
             
    plt.ylim((-350,0))
    plt.legend(names)
def save_currents(SP,US,VS,fn_out,nemo_t_out):
    
    coords = {
        "latitude": (("y_dim", "x_dim"), nemo_t_out.dataset.latitude.values),
        "longitude": (("y_dim", "x_dim"), nemo_t_out.dataset.longitude.values),
    }
    dims = ["y_dim", "x_dim"]
    nemo_t_out.dataset['Speed']=xr.DataArray(SP,coords=coords, dims=dims)
    nemo_t_out.dataset['U_unitvector']=xr.DataArray(US,coords=coords, dims=dims)
    nemo_t_out.dataset['V_unitvector']=xr.DataArray(VS,coords=coords, dims=dims)            

    nemo_t_out.dataset.to_netcdf(fn_out)
    
    
def regrid_currents(inputgrid,outputgrid):

    import xesmf as xe
    xesmf_ready = coast.xesmf_convert(inputgrid, outputgrid, output_grid_type = 'curvilinear',input_grid_type = 'curvilinear')
    regridder = xe.Regridder(xesmf_ready.input_grid,
                 xesmf_ready.output_grid, "bilinear")
    regridded_dataset = regridder(xesmf_ready.input_data)
    return regridded_dataset 
            
######################################################################
if __name__ == '__main__':
    fn_config_t_grid='./example_nemo_grid_t.json'
    fn_config_f_grid='./example_nemo_grid_f.json'
    fn_config_u_grid='./example_nemo_grid_u.json'
    fn_config_v_grid='./example_nemo_grid_v.json'


    #names,dpaths,DOMS,_  = coast. experiments(experiments='experiments_james.json')
    names=["EXP_MES_WAV_NTM"]
    dpaths=["/gws/nopw/j04/class_vol2/senemo/jdha/SHORT_TESTS/EXP_MES_WAV_NTM"]
    DOMS=["/gws/nopw/j04/class_vol2/senemo/jdha/SHORT_TESTS/EXP_MES_WAV_NTM/domain_cfg_r015-r010_007_004v2.nc"]       
    Unmean={}
    Untm={}
    Zmean={}
    Un={}
    Z={}
    ystart=1979
    ystop=1981
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
        fn_nemo_dat_u= coast.nemo_filename_maker(dpaths[i],ystart,ystop,grid='U',runtype='SENEMO') 
        fn_nemo_dat_v= coast.nemo_filename_maker(dpaths[i],ystart,ystop,grid='V',runtype='SENEMO')
        fn_nemo_dat_v=    fn_nemo_dat_v[iseason]   
        fn_nemo_dom=DOMS[i]
    #    fn_nemo_dat_u=dpaths[i]+'SENEMO_1m_19800101_19801231_grid_U_198007-198007.nc'
    #    fn_nemo_dat_v=dpaths[i]+'SENEMO_1m_19800101_19801231_grid_V_198007-198007.nc'
    
        recalc=True
        cmems=False            
        nemo_t = coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_t_grid)#,calc_bathy=True)
        if recalc:        
    
            nemo_f = coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_f_grid)#,calc_bathy=True)
            nemo_u = coast.Gridded(fn_data=fn_nemo_dat_u, fn_domain=fn_nemo_dom, config=fn_config_u_grid,multiple=True)
            nemo_v = coast.Gridded(fn_data=fn_nemo_dat_v, fn_domain=fn_nemo_dom, config=fn_config_v_grid,multiple=True)
        for ilme in lmelist:    
            LMENAM=A['DOMNAM'][ilme]

            print(LMENAM)
            x_min=A['x_min'][ilme]
            x_max=A['x_max'][ilme]
            y_min=A['y_min'][ilme]
            y_max=A['y_max'][ilme]        
                    
            j,i,_=nemo_t.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
            imin=min(i)
            imax=max(i)
            jmin=min(j)
            jmax=max(j)   
            jmin,jmax=[860,1015]
            imin,imax=[1080,1180]
            REGION=LMENAM
            REGION='NWS'
        
            if recalc:  
                nemo_f1=nemo_f.subset_as_copy(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
                nemo_u1=nemo_u.subset_as_copy(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
                nemo_v1=nemo_v.subset_as_copy(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
                nemo_t1=nemo_t.subset_as_copy(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
                mask=nemo_t1.dataset.bathymetry != 0
                Name=name+' '+SEASON+' '+YEARS+' '+REGION 
                SP,US,VS=mean_surface_circulation(nemo_u1,nemo_v1,nemo_t1,mask)
                plot_surface_circulation(SP,US,VS,nemo_t1,mask,Name)
                plt.savefig('../Figures/Circulation/Surface_Currents_' + Name.replace(' ','_')+'_v1.png')
                fn_out=("/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/Circulation/Surface_Currents_{0}.nc".format(Name)).replace(' ','_')
                nemo_t_out=coast.Gridded(fn_domain=fn_nemo_dom, config=fn_config_t_grid)
                nemo_t_out.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))
                save_currents(SP,US,VS,fn_out,nemo_t_out)

#
            if cmems:
                name="CMEMS_ORCA12"
                Name=name+' '+SEASON+' '+YEARS+' '+REGION   
                fn_cmems=("/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/Circulation/Surface_Currents_{0}.nc".format(Name)).replace(' ','_')
                cmems=coast.Gridded(fn_data=fn_cmems,config="")
                name=names[0]
                Name=name+' '+SEASON+' '+YEARS+' '+REGION
                fn_orca025=("/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/Circulation/Surface_Currents_{0}.nc".format(Name)).replace(' ','_')
                ORCA025=coast.Gridded(fn_data=fn_orca025,config="")
                cmems_on_ORCA025=regrid_currents(cmems,ORCA025)
            
    #    Un[name],Unmean[name],Z[name],Zmean[name]=flx_contour(nemo_f,nemo_u,nemo_v)
    
    
    
    
    
    #plot_flx_contour(Unmean,Zmean,names)
        
