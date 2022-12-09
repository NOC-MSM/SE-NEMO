import numpy as np
import matplotlib.pylab as plt

import pandas as pd
import sys
sys.path.insert(0, 'C:\\Users\\Jason.Goliath\\Documents\\GitHub\\COAsT\\')
import coast

LME_Data=np.load('../Data/LME_gridinfo_V4.npz')

LME_Clusters='../Data/LME_Clusters.csv'
clusters = pd.read_csv(LME_Clusters)
clusters = clusters.to_numpy()
cluster_names=clusters[:,0]
LME_numbers=clusters[:,1:]
Clusters={}
for iname,name in enumerate(cluster_names):
    lmes=LME_numbers[iname,:].astype('float')
    lmes=lmes[np.isfinite(lmes)].astype(int)
    Clusters[name]={}
    Clusters[name]['LMEs']=lmes

    Clusters[name]['limits']=np.array([np.min(LME_Data['x_min'][lmes-1]),
                                     np.max(LME_Data['x_max'][lmes-1]),
                                     np.min(LME_Data['y_min'][lmes-1]),
                                     np.max(LME_Data['y_max'][lmes-1])])



bathyname='../Data/eORCA025_bathy_meter.nc'
config='example_nemo_grid_t.json'
bathy=coast.Gridded(fn_data=bathyname,config=config)
D=np.ma.masked_where(bathy.dataset.Bathymetry.values==0,bathy.dataset.Bathymetry.values)
nx=D.shape[1]
plt.pcolormesh(D)
for name in cluster_names:
    lims=Clusters[name]['limits']
    j, i, _ = bathy.find_j_i_list(lon=[lims[0], lims[1], lims[1], lims[0]], lat=[lims[2], lims[2],lims[3], lims[3] ])

    imin = min(i)
    imax = max(i)
    jmin = min(j)
    jmax = max(j)
    if i[0] < i[1]:
      plt.plot([imin,imax,imax,imin,imin],[jmin,jmin,jmax,jmax,jmin])
    else:
      print(name)
      plt.plot([imax, nx, np.nan, nx, imax, imax,np.nan,
                0,imin,imin,0],
               [jmin, jmin,np.nan, jmax, jmax, jmin,np.nan,
                jmin,jmin,jmax,jmax])
      

