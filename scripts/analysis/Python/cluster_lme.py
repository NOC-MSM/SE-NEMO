import numpy as np
import matplotlib.pylab as plt

import pandas as pd
import sys
sys.path.insert(0, 'C:\\Users\\Jason.Goliath\\Documents\\GitHub\\COAsT\\')
sys.path.insert(0,'/home/users/jholt/Git/COAsT/')

import coast
import scipy.io
LME_Data=np.load('../Data/LME_gridinfo_V4.npz')
a=scipy.io.loadmat('../Data/ORCA025_ROAM_GLB_LMEmaskV4.mat')
LME_mask=a['LME_mask'][:,:].T
LME_Clusters='../Data/LME_Clusters.csv'
LME_Clusters_out='../Data/LME_Clusters_eORCA025.csv'

clusters = pd.read_csv(LME_Clusters)
clusters = clusters.to_numpy()
cluster_names=clusters[:,0]
LME_numbers=clusters[:,1:]
Clusters={}

bathyname='../Data/eORCA025_bathy_meter.nc'
bathyname='/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/config/domain_cfg.nc'

config='example_nemo_grid_t.json'
bathy=coast.Gridded(fn_data=bathyname,config=config)
D=np.ma.masked_where(bathy.dataset.bathymetry.values==0,bathy.dataset.bathymetry.values).squeeze()
lon=bathy.dataset.longitude
lat=bathy.dataset.latitude


nx=D.shape[1]
J_offset=186 #account for extra rows in eORCA if data is made for normal ORCA

for iname,name in enumerate(cluster_names):
    lmes=LME_numbers[iname,:].astype('float')
    lmes=lmes[np.isfinite(lmes)].astype(int)
    Clusters[name]={}
    Clusters[name]['LMEs']=lmes
    lims=[99999,-99999,99999,-99999]
    X=np.array([])
    Y=np.array([])

    for lme in lmes:
        a=np.where(LME_mask==lme)
        X=np.append(X,a[1])
        Y=np.append(Y,a[0])

 


    lims=np.array([np.min(X),np.max(X),np.min(Y),np.max(Y)]).astype(int)
    if name == 'S Asia':
        lims[0]=1280
        lims[1]=122
    if name == 'Barents Sea':
        lims[3]=1021
    if name == 'E Canadian Arctic':
        lims[0]=784
    if name == 'W Canadian Arctic':
        lims[1]=700
    if name == 'Eurasian Arctic':
        lims[1]=341
    if name == 'Nordic Seas':
        lims[0] = 981        
    Clusters[name]['limits']=lims+[0,0,J_offset,J_offset]
    
    DD=coast.Gridded(fn_data=bathyname,config=config,lims=lims+[0,0,J_offset,J_offset])
    x=DD.dataset.longitude.values
    if np.max(x)-np.min(x)>350:
        x[x<0]=x[x<0]+360
    y=DD.dataset.latitude.values
    x_masked=np.ma.masked_where(DD.dataset.bottom_level.values[0,:,:]==0,x)
    y_masked=np.ma.masked_where(DD.dataset.bottom_level.values[0,:,:]==0,y)    
    xylims=np.array([np.min(x_masked),np.max(x_masked),np.min(y_masked),np.max(y_masked)])
    Clusters[name]['xylimits']=xylims  
nname=iname+1
  
plt.pcolormesh(D[J_offset:,:])
plt.ylabel('j-189')
Lims=np.zeros((nname,4))
xyLims=np.zeros((nname,4))
for iname,name in enumerate(cluster_names):
    lims=Clusters[name]['limits']
    xylims= Clusters[name]['xylimits']
    imin = lims[0]
    imax = lims[1]
    jmin = lims[2]-J_offset
    jmax = lims[3]-J_offset
    if name !='S Asia':
      plt.plot([imin,imax,imax,imin,imin],[jmin,jmin,jmax,jmax,jmin])
    else:
      print(name)
      plt.plot([imin, nx, np.nan, nx, imin, imin,np.nan,
                0,imax,imax,0],
               [jmin, jmin,np.nan, jmax, jmax, jmin,np.nan,
                jmin,jmin,jmax,jmax])
    Lims[iname,:]=lims
    xyLims[iname,:]=xylims
      
LME_Data['DOMNAM'][63]=LME_Data['DOMNAM'][53]      
PD={}
PD['Cluster']=cluster_names
PD['i min']=Lims[:,0]
PD['i max']=Lims[:,1]
PD['j min']=Lims[:,2]
PD['j max']=Lims[:,3]
PD['x min']=xyLims[:,0]
PD['x max']=xyLims[:,1]
PD['y min']=xyLims[:,2]
PD['y max']=xyLims[:,3]

for ilme in range(LME_numbers.shape[1]  ):
    NN=[]
    for iname,name in enumerate(cluster_names):     
        try:
            NN.append(LME_Data['DOMNAM'][int(LME_numbers[iname,ilme])-1])     
        except:
            NN.append('')
    PD['LME {0}'.format(ilme+1)]=NN

df=pd.DataFrame(PD)
df.to_csv(LME_Clusters_out)    





