import numpy as np
import matplotlib.pylab as plt
import pickle
# %%
from matplotlib import cm
cmap0=cm.get_cmap('BrBG_r',lut=32)
plt.ion
LME_gridinfo=np.load('../Data/LME_gridinfo_V4.npz')
nlme=66
ind_FC_mean={}
ind_FC_max={}
ind_FC_min={}

keys=['sst','sss','ssso','pea','nitrate','O2_bot','dsss']
for key in keys:
   ind_FC_max[key] = np.zeros(nlme)
   ind_FC_mean[key] = np.zeros(nlme)
   ind_FC_min[key]  = np.zeros(nlme)


domain_outpath = '/home/users/jholt/work/SENEMO/SENEMO_FC/'
ind_sc = {}
ind_mean = {}
ind_max = {}
ind_min = {}

nsl = 20

for ilme in range(66):#np.concatenate([np.arange(0,16),[23]]):
   #%%
   lme_name = LME_gridinfo['DOMNAM'][ilme]
   inname = f"{domain_outpath}physical_indicators_bgc_{lme_name}.p"
   try:
   #%%
      with open(inname,'rb' ) as f:
         indicators=pickle.load(f)
      indicators['dsss'] =  indicators['sss'] / (indicators['ssso'] - indicators['sss'])

      it0 = np.arange(0,nsl)
      it1 = it0 + (2070-nsl+1-1980)


      for key in list(indicators.keys()):
         var = indicators[key]
         nyears = int((var.shape[0]/12))
         ind_sc[key]=indicators[key].reshape(nyears,12)
         ind_mean[key,ilme]=np.mean(ind_sc[key],axis=1)
         ind_max[key,ilme] = np.max(ind_sc[key], axis=1)
         ind_min[key, ilme] = np.min(ind_sc[key], axis=1)

         ind_FC_mean[key][ilme] = (np.mean(ind_mean[key,ilme][it1]) - np.mean(ind_mean[key,ilme][it0]))/ np.std(ind_mean[key,ilme][it0])
         ind_FC_max[key][ilme] =  (np.mean(ind_max[key,ilme][it1]) - np.mean(ind_max[key,ilme][it0])) / np.std(ind_max[key,ilme][it0])
         ind_FC_min[key][ilme] = (np.mean(ind_min[key, ilme][it1]) - np.mean(ind_min[key, ilme][it0])) / np.std(ind_min[key, ilme][it0])

      ind_FC_mean['sst'][ilme] = (np.mean(ind_mean['sst',ilme][it1]) - np.mean(ind_mean['sst',ilme][it0]))
      ind_FC_max['sst'][ilme] =  (np.mean(ind_max['sst',ilme][it1]) - np.mean(ind_max['sst',ilme][it0]))
   #%%
   except:
      print('missing',lme_name)

#%%
LME_gridinfo=np.load('../Data/LME_gridinfo_V4.npz')
names = LME_gridinfo['DOMNAM']

nind=len(list(indicators.keys()))
nlme=66
ind_grid=np.zeros((nlme,5))
ind_grid[:,0] = ind_FC_min['O2_bot']
ind_grid[:,1] = ind_FC_max['nitrate']
ind_grid[:,2] = ind_FC_mean['sst']
ind_grid[:,3] = ind_FC_mean['dsss']
ind_grid[:,4] = ind_FC_max['pea']
labs = ["""O2/$\sigma$
Bot Min""","""Nit/$\sigma$
Surf Max""","""SST oC
Mean""","""DSSS/$\sigma$
Mean""","""PEA/$\sigma$
Max"""]
plt.ion()

plt.figure(figsize=(8.27,11.69))
plt.pcolormesh(ind_grid,vmin=-4,vmax=4,cmap=cmap0)
ax=plt.gca()
ax.set_yticks(np.arange(0.5,nlme+0.5),names[:nlme],fontsize=8)
ax.set_xticks(np.arange(0.5,len(labs)+0.5),labs,fontsize=8)

ax.set_position([0.4,.15,0.5,0.8])
ax.set_title("""Change 2051-2070 - 1980-1999
Average over Coastal LMEs <500m""")
#plt.colorbar(orientation='horizontal')

ax_cb = plt.axes(position = [.4,.08,0.5,0.02])
plt.colorbar(ax=ax,cax=ax_cb, orientation = 'horizontal')
plt.savefig('../Figures/Inds_LME_cnrm.png',dpi=300)


#%%
key='nitrate'
fig,axs = plt.subplots(nrows=8,ncols=9,figsize=[11.69,8.27])
axs=axs.ravel()
t=np.arange(1980,2071)
for ilme in np.arange(72):

   if (key,ilme) in list(ind_mean.keys()):
      lme_name = LME_gridinfo['DOMNAM'][ilme]
      axs[ilme].plot(t,ind_mean[key,ilme][:])
      axs[ilme].set_title(f"LME: {ilme+1}",fontsize=8)
      axs[ilme].tick_params(axis='both',labelsize=8)
      if ilme <= 56:
         axs[ilme].set_xticks([])
      else:
         axs[ilme].set_xticks([2000,2050],labels=['2000','2050'],fontsize=8)
   else:
      print('missing')
      axs[ilme].remove()
