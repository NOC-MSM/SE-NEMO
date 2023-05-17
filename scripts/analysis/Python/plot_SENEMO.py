
import sys
sys.path.insert(0,'/login/jholt/work/Git/COAsT/')
import matplotlib.pylab as plt
import coast
import numpy as np
#cmap0=cm.get_cmap('BrBG_r',lut=16)
cmap0.set_bad(color=[0.65,0.65,0.65])
#%%
fn_nemo_dom='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/domcfg_eORCA025_v2.nc'
fn_nemo_dat1='/work/jholt/JASMIN//SENEMO/NOTIDE/SENEMO_1m_19800101_19801231_grid_T_198001-198001.nc'
fn_nemo_dat2=  '/work/jholt/JASMIN//SENEMO/TIDE/SENEMO_1m_19800101_19801231_grid_T_198001-198001.nc'

fn_nemo_dat1='/work/jholt/JASMIN//SENEMO/JDHA/EXP_ZPS/SENEMO_1M/SENEMO_1m_19800101_19801231_grid_T_198008-198008.nc'
fn_nemo_dat2='/work/jholt/JASMIN//SENEMO/JDHA/EXP_SZT39_TAPER/SENEMO_1M/SENEMO_1m_19800101_19801231_grid_T_198008-198008.nc'
fn_nemo_dom='/work/jholt/JASMIN//SENEMO/JDHA/EXP_SZT39_TAPER/domain_cfg_ztaper_match.nc'


fn_nemo_dat_w1='/work/jholt/JASMIN//SENEMO/NOTIDE/SENEMO_1m_19800101_19801231_grid_W_198001-198001.nc'
fn_nemo_dat_w2=  '/work/jholt/JASMIN//SENEMO/TIDE/SENEMO_1m_19800101_19801231_grid_W_198001-198001.nc'

fn_config_t_grid='/login/jholt/work/Git/COAsT/config/example_nemo_grid_t.json'
fn_config_w_grid='/login/jholt/work/Git/COAsT/config/example_nemo_grid_w.json'


nemo_t1 = coast.Gridded(fn_data = fn_nemo_dat1, fn_domain = fn_nemo_dom, config=fn_config_t_grid)
nemo_t2 = coast.Gridded(fn_data = fn_nemo_dat2, fn_domain = fn_nemo_dom, config=fn_config_t_grid)


nemo_w1 = coast.Gridded(fn_data = fn_nemo_dat_w1, fn_domain = fn_nemo_dom, config=fn_config_w_grid)
nemo_w2 = coast.Gridded(fn_data = fn_nemo_dat_w2, fn_domain = fn_nemo_dom, config=fn_config_w_grid)



SSS1=nemo_t1.dataset.variables['so_abs'].values[0,0,:,:]
SSS2=nemo_t2.dataset.variables['so_abs'].values[0,0,:,:]
SST1=nemo_t1.dataset.variables['thetao_con'].values[0,0,:,:]
SST2=nemo_t2.dataset.variables['thetao_con'].values[0,0,:,:]


plt.pcolormesh(SSS2)

#%%
j=660
i=200


j=950
i=1140

i=565
j=685

avh1=nemo_w1.dataset.variables['difvho'].values[0,:,j,i]
avh2=nemo_w2.dataset.variables['difvho'].values[0,:,j,i]
T1=nemo_t1.dataset.variables['thetao_con'].values[0,:,j,i]
T2=nemo_t2.dataset.variables['thetao_con'].values[0,:,j,i]

Z1=nemo_t1.dataset.coords['depth_0'].values[:,j,i]



plt.subplot(1,2,1)
plt.plot(T1,-Z,T2,-Z,'--')
plt.legend(['No Tide tmp','Tide tmp'])
# =============================================================================
# plt.title('Profiles at j,i={0},{1}'.format(j,i))
# =============================================================================
plt.subplot(1,2,2)
plt.plot(np.log10(avh1),-Z,np.log10(avh2),-Z,'--')
plt.legend(['No Tide avt (log10)','Tide avt (log10)'])

plt.savefig('/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT//ORCA025-SE-NEMO/prof_j{0}_i{1}_1980_tide-notide.png'.format(j,i))

#%%%
plt.close('all')
dpath='/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT/ORCA025-SE-NEMO/'


EXPNAMS=['EXP_MES_TIDE','EXP_MES_NOTIDE','TIDE','NOTIDE']
#EXPNAMS2=['TIDE','NOTIDE']
pea_ns={}
nemo={}

for EXPNAM in EXPNAMS:
 fname=dpath+'ORCA025-SE-NEMO_1980_1985_' +EXPNAM+ '_SST_SSS_PEA_MonClimate.nc'   
 nemo[EXPNAM]=coast.Gridded(fn_data=fname)
#for EXPNAM in EXPNAMS2:
# fname=dpath+'ORCA025-SE-NEMO/ORCA025-SE-NEMO__TSclim_'+ EXPNAM +'1980_2011.nc'   
# nemo[EXPNAM]=coast.Gridded(fn_data=fname)
SST1=nemo[EXPNAMS[0]].dataset.variables['SSTy'][0,860:1000,1080:1180]
for EXPNAM in EXPNAMS:
 pea=np.mean(nemo[EXPNAM].dataset.variables['PEAy'].values[5:9,860:1000,1080:1180],axis=0)
     
 pea_ns[EXPNAM]=np.ma.masked_where(np.isnan(SST1),pea)
#for EXPNAM in EXPNAMS2:
# pea=nemo[EXPNAM].dataset.variables['PEAy'].values[1080:1180,860:1000,7].T
# pea_ns[EXPNAM]=np.ma.masked_where(np.isnan(SST1),pea)
for EXPNAM in EXPNAMS: #+EXPNAMS2:

 plt.figure()
 plt.pcolormesh(pea_ns[EXPNAM],vmin=0,vmax=200,cmap=cmap0)
 plt.colorbar(orientation='vertical')
 plt.title('PEA ' +EXPNAM)
 plt.savefig('/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT//ORCA025-SE-NEMO/PEA_NWS_{0}.png'.format(EXPNAM))





