import sys
import socket
isliv = 'livljobs' in socket.gethostname()
if isliv:
 sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
else:
 sys.path.insert(0,'/home/users/jholt/Git/COAsT/')

import matplotlib.pylab as plt
import coast
import numpy as np
import pickle
import pandas as pd
outname='/home/users/jholt/work/SENEMO/ASSESSMENT/ORCA025-SE-NEMO/SL_std_SENEMO.p'


names,dpaths,DOMS,_  = coast.experiments(experiments='../Python/experiments.json')
 
ystart=1990
ystop=2019


fn_config_t_grid='../Config/senemo_grid_t.json'
LME_Clusters='../Data/LME_Clusters_eORCA025.csv'
clusters = pd.read_csv(LME_Clusters)

#%%
ssh_std={}
ssh_mn={}
slmean=np.zeros((clusters.values.shape[0],3))
for iexp in [0,1,2]:#range(3):
    fn_nemo_dat= coast.nemo_filename_maker(dpaths[iexp],ystart,ystop)
    for icluster in range(clusters.values.shape[0]):
      #if icluster != 10:  #COAST subsettign across long wrap doesnt work here??
        print(iexp,icluster)
        lims=np.array(clusters.values[icluster,2:6],dtype=int)    
        nemo = coast.Gridded(fn_data= fn_nemo_dat, fn_domain= DOMS[iexp],config=fn_config_t_grid,multiple=True,lims=lims)
        ssh_std[icluster,iexp]=np.std(nemo.dataset.ssh.values[:,:,:],axis=0)
        nt = nemo.dataset.dims['t_dim']
        DX=np.repeat(nemo.dataset.e1.values[np.newaxis,:,:],nt,axis=0)
        DY=np.repeat(nemo.dataset.e1.values[np.newaxis,:,:],nt,axis=0)
        
        slmean[icluster,iexp]=np.nansum(nemo.dataset.ssh.values[:,:,:]*DX*DY)/np.nansum(DX*DY)
        ssh_mn[icluster,iexp]=np.mean(nemo.dataset.ssh.values[:,:,:],axis=0) - slmean[icluster,iexp]

#%%
with open(outname,'wb' ) as f:
   A={}
   A['ssh_std']=ssh_std
   A['ssh_mn']=ssh_mn
   A['slmean']=slmean
   pickle.dump(A,f)

#np.save(outname+'_std1.npy',ssh_std)
        #,ssh_std=ssh_std,ssh_mn=ssh_mn,slmean=slmean)
#plt.pcolormesh(ssh_std[0,0])
#plt.colorbar(orientation='vertical')
#%%

with open(outname,'rb' ) as f:
    A=pickle.load(f)
ssh_std=A['ssh_std']
ssh_mn=A['ssh_mn']
slmean=A['slmean']
iexp=0
import plot_surfacefield_cluster as pl_clst
x={}
y={}
dssh_std={}
dssh_mn={}
dssh_gmn={}


for icluster in range(clusters.values.shape[0]):
      #if icluster != 10:  #COAST subsettign across long wrap doesnt work here??
        print(iexp,icluster)
        lims=np.array(clusters.values[icluster,2:6],dtype=int)    
        nemo = coast.Gridded(fn_domain= DOMS[iexp],config=fn_config_t_grid,lims=lims)
        x[icluster]=nemo.dataset.longitude.values
        y[icluster]=nemo.dataset.latitude.values
        dssh_std[icluster]=ssh_std[icluster,0]-ssh_std[icluster,2]   
        dssh_mn[icluster]=ssh_mn[icluster,0]-ssh_mn[icluster,2]   
        dssh_gmn[icluster]=ssh_mn[icluster,0]-ssh_mn[icluster,2]+slmean[icluster,0]-slmean[icluster,2] 

#%%
vmin=-.1
vmax=.1
Title='Change in SSH std (m)'
Figname='../Figures/SSH_std_SENEMO-ZPS.png'
pl_clst.cluster_plot(x,y,dssh_std,vmin,vmax,Title,Figname)

#%%
vmin=-.2
vmax=.2
Title='Change in SSH local mean (m)'
Figname='../Figures/SSH_lmean_SENEMO-ZPS.png'
pl_clst.cluster_plot(x,y,dssh_mn,vmin,vmax,Title,Figname)
#%%
vmin=-.2
vmax=.2
Title='Change in SSH mean (m)'
Figname='../Figures/SSH_gmean_SENEMO-ZPS.png'
pl_clst.cluster_plot(x,y,dssh_gmn,vmin,vmax,Title,Figname)



