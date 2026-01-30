import numpy as np
import matplotlib.pylab as plt
import pickle
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
# %%

def lightcolormap(Np, nt):
   cmap0 = cm.get_cmap('BrBG_r', lut=Np + nt * 2)
   colors = cmap0(np.arange(cmap0.N))
   colors1 = colors[nt:cmap0.N - nt]
   cmap1 = LinearSegmentedColormap.from_list('cmap1', colors1, cmap0.N - nt * 2)
   return cmap1
cmap0=cm.get_cmap('BrBG_r',lut=32)
cmap0=lightcolormap(32,2)
plt.ion
LME_gridinfo=np.load('../Data/LME_gridinfo_V4.npz')
nlme=66
ind_FC_mean={}
ind_FC_max={}
ind_FC_min={}

keys=['sst','nbt','sss','ssso','pea','nitrate','O2_bot','dsss','Qu','Qv','Q','nitrate_o','netpp','mld','mldo','ice','runoff']
iwant_lme=np.array([5,6,7,8,9,12,14,15,16,17,18,19,20,21,22,24,25,26,27,28,29,59,60,63])-1
iwant_lme=np.array([14,15,16,17,12,5,6,7,8,9,63,18,19,59,20,21,60,24,22,25,26,27,28,29])-1

nlme=len(iwant_lme)
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
ESM='CNRM'
SSP='ssp370'

for icount,ilme in enumerate(iwant_lme):#(66):#np.concatenate([np.arange(0,16),[23]]):
   #%%
   lme_name = LME_gridinfo['DOMNAM'][ilme]
   inname     = f"{domain_outpath}physical_indicators_{ESM}_{SSP}_{lme_name}_200m.p"
   inname_npp = f"{domain_outpath}physical_indicators_{ESM}_{SSP}_{lme_name}_200m.p"
   #if ilme==21:
   #   inname = f"{domain_outpath}physical_indicators_{ESM}_{SSP}_{lme_name}_200m.p"
   #   inname_npp = f"{domain_outpath}physical_indicators_{ESM}_{SSP}_{lme_name}_200m.p"
   #try:
   if True:
   #%%
      with open(inname,'rb' ) as f:
         indicators=pickle.load(f)
      with open(inname_npp,'rb' ) as f:
         indicators_npp=pickle.load(f)
      indicators['dsss'] = indicators['ssso'] - indicators['sss']
      indicators['Q'] =  indicators['Qu'] + indicators['Qv']
      indicators['netpp']=indicators_npp['netpp']*365/1000
      it0 = np.arange(0,nsl)
      it1 = it0 + (2070-nsl+1-1980)


      for key in list(indicators.keys()):
         var = indicators[key]
         nyears = int((var.shape[0]/12))
         ind_sc[key]=indicators[key].reshape(nyears,12)
         ind_mean[key,ilme]=np.mean(ind_sc[key],axis=1)
         ind_max[key,ilme] =np.max(ind_sc[key], axis=1)
         ind_min[key, ilme] = np.min(ind_sc[key], axis=1)
         if 'runoff' in key:
            ind_mean[key, ilme] = ind_sc[key][:,0]
         ind_FC_mean[key][icount] = (np.mean(ind_mean[key,ilme][it1]) - np.mean(ind_mean[key,ilme][it0]))/ np.std(ind_mean[key,ilme][it0])
         ind_FC_max[key][icount] =  (np.mean(ind_max[key,ilme][it1]) - np.mean(ind_max[key,ilme][it0])) / np.std(ind_max[key,ilme][it0])
         ind_FC_min[key][icount] = (np.mean(ind_min[key, ilme][it1]) - np.mean(ind_min[key, ilme][it0])) / np.std(ind_min[key, ilme][it0])

      ind_FC_mean['sst'][icount] = (np.mean(ind_mean['sst',ilme][it1]) - np.mean(ind_mean['sst',ilme][it0]))
      ind_FC_max['sst'][icount] =  (np.mean(ind_max['sst',ilme][it1]) - np.mean(ind_max['sst',ilme][it0]))
      ind_FC_mean['nbt'][icount] = (np.mean(ind_mean['nbt', ilme][it1]) - np.mean(ind_mean['nbt', ilme][it0]))
      ind_FC_max['nbt'][icount] = (np.mean(ind_max['nbt', ilme][it1]) - np.mean(ind_max['nbt', ilme][it0]))
   #%%
   else:
   #except:
      print('missing',lme_name)

#%%
LME_gridinfo=np.load('../Data/LME_gridinfo_V4.npz')
names=[]
for icount in iwant_lme:
   names.append(LME_gridinfo['DOMNAM'][icount].replace('_',' '))

nind=len(list(indicators.keys()))

ind_grid=np.zeros((nlme,14))
ind_grid[:,0] = ind_FC_mean['netpp']
ind_grid[:,1] = ind_FC_mean['O2_bot']

ind_grid[:,2] = ind_FC_max['nitrate']
ind_grid[:,3] = ind_FC_max['nitrate_o']
ind_grid[:,4] = ind_FC_mean['Q']
ind_grid[:,5] = ind_FC_mean['sst']
ind_grid[:,6] = ind_FC_mean['nbt']
ind_grid[:,7] = ind_FC_mean['dsss']
ind_grid[:,8] = ind_FC_max['pea']
ind_grid[:,9] = ind_FC_max['mld']
ind_grid[:,10] = ind_FC_max['mldo']
ind_grid[:,11] = ind_FC_min['mld']
ind_grid[:,12] = ind_FC_mean['runoff']
ind_grid[:,13] = ind_FC_min['ice']




labs = ["""netPP/$\sigma$
Mean""","""O2/$\sigma$
Bot Min""","""Nit/$\sigma$
Surf Max""",
"""Nit O/$\sigma$
Surf Max""",
"""Q/$\sigma$
Mean""",
"""SST oC
Mean""",
"""NBT oC
Mean""",
"""$\Delta$S/$\sigma$
Mean""",
"""PEA/$\sigma$
Max""","""MLD/$\sigma$
Max""","""MLD O/$\sigma$
Max""","""MLD/$\sigma$
Min""",
"""Runoff/$\sigma$
Mean""",
"""Ice/$\sigma$
Min"""
        ]
plt.ion()

#plt.figure(figsize=(8.27,11.69))
plt.figure(figsize=(11.69,8.27))

plt.pcolormesh(ind_grid,vmin=-4,vmax=4,cmap=cmap0)
ax=plt.gca()
ax.set_position([0.325,.15,0.65,0.8])
ax.set_yticks(np.arange(0.5,nlme+0.5),names[:nlme],fontsize=8)
ax.set_xticks(np.arange(0.5,len(labs)+0.5),labs,fontsize=8)


ax.set_title(f"""Change 2051-2070 - 1980-1999 {ESM}
Average over Coastal LMEs <200m""")
#plt.colorbar(orientation='horizontal')

ax_cb = plt.axes(position = [.4,.08,0.5,0.02])
plt.colorbar(ax=ax,cax=ax_cb, orientation = 'horizontal')
plt.savefig(f'../Figures/Inds_LME_{ESM}_200m.png',dpi=300)


#%%
key='netpp'
ylabel='netpp (gCm$^{-2}$yr$^{-1}$)'
#key=('O2_bot')
#ylabel='O$_{2}$ (mmol m$^{-3}$)'
fig,axs = plt.subplots(figsize=[11.69,8.27],nrows=5,ncols=5)#,nrows=8,ncols=9,)
axs=axs.ravel()
t=np.arange(1980,2071)
for icount,ilme in enumerate(iwant_lme):#np.arange(72):
   lme_name = LME_gridinfo['DOMNAM'][ilme].replace('_',' ')
   if (key,ilme) in list(ind_mean.keys()):

      axs[icount].plot(t,ind_mean[key,ilme][:])
      axs[icount].set_title(f"LME: {ilme+1}",fontsize=8)
      axs[icount].set_title(lme_name, fontsize=8)

      axs[icount].tick_params(axis='both',labelsize=8)
      axs[icount].grid()
      if icount <= 18:
         axs[icount].set_xticks([])
         axs[icount].set_xticks([2000,2020,2040,2060],labels=[],fontsize=8)


      else:
         axs[icount].set_xticks([2000,2020,2040,2060],labels=['2000','2020','2040','2060'],fontsize=8)
      if sum(icount == np.array([0,5,10,15,20])):
         axs[icount].set_ylabel(ylabel,fontsize=8)
   else:
      print('missing',key,ilme,lme_name)
      axs[icount].remove()
axs[icount+1].remove()
plt.savefig(f'../Figures/Series_LME_{key}_{ESM}_200m.png',dpi=300)