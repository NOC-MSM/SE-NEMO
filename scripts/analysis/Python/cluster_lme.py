import numpy as np
import matplotlib.pylab as plt

import pandas as pd
import sys
sys.path.insert(0, 'C:\\Users\\Jason.Goliath\\Documents\\GitHub\\COAsT\\')
import coast
import scipy.io
LME_Data=np.load('../Data/LME_gridinfo_V4.npz')
a=scipy.io.loadmat('../Data/ORCA025_ROAM_GLB_LMEmaskV4.mat')
LME_mask=a['LME_mask'][:,:].T
LME_Clusters='../Data/LME_Clusters.csv'
clusters = pd.read_csv(LME_Clusters)
clusters = clusters.to_numpy()
cluster_names=clusters[:,0]
LME_numbers=clusters[:,1:]
Clusters={}

bathyname='../Data/eORCA025_bathy_meter.nc'
config='example_nemo_grid_t.json'
bathy=coast.Gridded(fn_data=bathyname,config=config)
D=np.ma.masked_where(bathy.dataset.Bathymetry.values==0,bathy.dataset.Bathymetry.values)
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
    for lme in lmes:
        a=np.where(LME_mask==lme) 

   
        lims=np.array([np.min([lims[0],np.min(a[1])]),
               np.max([lims[1],np.max(a[1])]),
               np.min([lims[2],np.min(a[0])]),
               np.max([lims[3],np.max(a[0])])    
             ])

    Clusters[name]['limits']=lims+[0,0,J_offset,J_offset]

plt.pcolormesh(D)
for name in cluster_names:
    lims=Clusters[name]['limits']

    imin = lims[0]
    imax = lims[1]
    jmin = lims[2]
    jmax = lims[3]
    if name !='S Asia':
      plt.plot([imin,imax,imax,imin,imin],[jmin,jmin,jmax,jmax,jmin])
    else:
      print(name)
      plt.plot([imax, nx, np.nan, nx, imax, imax,np.nan,
                0,imin,imin,0],
               [jmin, jmin,np.nan, jmax, jmax, jmin,np.nan,
                jmin,jmin,jmax,jmax])
      

