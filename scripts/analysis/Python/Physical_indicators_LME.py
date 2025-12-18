

import numpy as np
import scipy.io
import matplotlib.pylab as plt
import time
import xarray as xr
import pickle
import sys
args=sys.argv
if len(sys.argv) >= 2:
    iwant_lme = int(args[1])
else:
    iwant_lme = 0
print ('LME_list',iwant_lme)


sys.path.insert(0,'/home/users/jholt/Git/COAsT/')
#needs branch     feature/535_stratification_diag
import coast

try:
    plt.figure()
except:
    print('error in matplotlib')
####
def is_leap(year):
    """
    Checks if a given year is a leap year.
    Args:
        year: The year to check (integer).
    Returns:
        True if the year is a leap year, False otherwise.
    """
    # Rule 1: Must be divisible by 4
    if (year % 4) == 0:
        # Rule 2: If divisible by 100, import xarray as xrmust also be divisible by 400
        if (year % 100) == 0:
            if (year % 400) == 0:
                return True
            else:
                return False
        else:
            return True
    else:
        return False
####
LME_gridinfo=np.load('../Data/LME_gridinfo_V4.npz')
a=scipy.io.loadmat('../Data/ORCA025_ROAM_GLB_LMEmaskV4.mat')
LME_mask=a['LME_mask'][:,:].T
J_offset=186 #account for extra rows in eORCA if data is made for normal ORCA
names,dpaths,DOMS,_,year_start,year_stop  = coast.experiments(experiments='experiments_FC_Indicators.json')

grid = 'T'
#%%
vars = ['sst','sss','ssso','pea']
varname =['temperature','salinity','salinity','pea']

vars_bgc = ['nitrate','nitrate_o',O2_bot']
varname_bgc =['N3_n','N3_n','O2_bot']

#%%
ystart0 = 1980
ilme=21
nlme=66
ntmax=(max(year_stop)-min(year_start)+1)*12

#%%
lme_list = {}
for i in range (10):
   lme_list[i] = np.arange(6) + i*6

lme_list[11] = np.arange(61,65)
lme_list[0]=np.array([22,33,34,35,36,38,39,40,41,42,59,60,65,66])-1
#%%
#iwant_lme = 11
for ilme in lme_list[iwant_lme]: #range(11,61):#,23]:
    indicators = {}

    for var in vars:
        indicators[var] = np.zeros((ntmax))
    for var in vars_bgc:
        indicators[var] = np.zeros((ntmax))

    print(LME_gridinfo['DOMNAM'][ilme])
    #%%
    LME_name = LME_gridinfo['DOMNAM'][ilme]
    imin = LME_gridinfo['i_min'][ilme]-1
    imax = LME_gridinfo['i_max'][ilme]+1
    jmin = LME_gridinfo['j_min'][ilme] + J_offset-1
    jmax = LME_gridinfo['j_max'][ilme] + J_offset+1
    jmin0 = LME_gridinfo['j_min'][ilme]-1
    jmax0 = LME_gridinfo['j_max'][ilme]+1
    #%%
    nemos={}
    for i in [0,1,2,3]:#, 1, 2, 3]:
        EXPNAM = names[i]
        ystart = year_start[i]
        ystop = year_stop[i]
        for year in range(ystart,ystop+1):
            print(year,'opening',EXPNAM)
            t0 = time.time()
            domain_datapath = dpaths[i]
            # make list of filenames
            fn_nemo_dat = coast.nemo_filename_maker(domain_datapath, ystart, ystop)
            fn_nemo_dat = []
            if 'CNRM' in EXPNAM:
                ESM = 'CNRM'
            else:
                ESM = 'GFDL'
            if 'hist' in EXPNAM:
                SSP = 'hist'
            else:
                SSP = 'ssp370'

            if 'bgc' in EXPNAM:
                new_name = f"{domain_datapath}/SE_{ESM}_subBGC_{SSP}_{year}.nc"
                fn_nemo_dat.append(new_name)
            else:
                new_name = f"{domain_datapath}/SE_{ESM}_forTransp_T_{SSP}_{year}.nc"
                fn_nemo_dat.append(new_name)

                if False:
                    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
                    if is_leap(year):
                        days[1] = 29
                    for Month in range(1, 13):
                        day_stop = days[Month - 1]
                        new_name = f"{domain_datapath}/{EXPNAM}_1m_{year}{Month:02}01_{year}{Month:02}{day_stop:02}_grid_{grid}_{year}{Month:02}-{year}{Month:02}.nc"
                        fn_nemo_dat.append(new_name)

            # Provide a domain.cfg file
            fn_nemo_dom = DOMS[i]

            # Provide a config file
            fn_config_t_grid = '../Config/senemo_grid_t.json'

            # input datasets
            lims=[imin,imax+1,jmin,jmax+1]
            nemo = coast.Gridded(fn_data=fn_nemo_dat, fn_domain=fn_nemo_dom, config=fn_config_t_grid,
                                 multiple=True,no_depths=True,lims=lims)
            t1 = time.time()
            print(t1-t0)


            # %%

            ntimes = nemo.dataset.sizes['t_dim']

            mask = nemo.dataset.variables['bottom_level'].values != 0
            if len(mask.shape) == 2:
                mask = np.repeat(mask[np.newaxis, :, :], ntimes, axis=0)
            lme_mask = mask * LME_mask[jmin0:jmax0 + 1, imin:imax + 1] == ilme + 1
            #lme_mask = np.repeat(lme_mask[np.newaxis, :, :], ntimes, axis=0)

            Depth_lim = 500
            Depth=nemo.dataset['bathymetry'].values
            Dmask = Depth <= Depth_lim
            Dmasko = Depth > Depth_lim
            #DX=np.repeat(nemo.dataset['e1'].values[np.newaxis, :, :], ntimes, axis=0)
            #DY=np.repeat(nemo.dataset['e2'].values[np.newaxis, :, :], ntimes, axis=0)
            DX=nemo.dataset['e1']
            DY=nemo.dataset['e2']
            if len(Dmask.shape)==2:
                Dmask = np.repeat(Dmask[np.newaxis, :, :], ntimes, axis=0)
                Dmasko = np.repeat(Dmasko[np.newaxis, :, :], ntimes, axis=0)
                DX = np.repeat(DX.values[np.newaxis, :, :], ntimes, axis=0)
                DY = np.repeat(DY.values[np.newaxis, :, :], ntimes, axis=0)

            area = np.sum(np.sum(Dmask*lme_mask*DX*DY,axis=2),axis=1)
            #%%
#Define external Points
            clme_mask = Dmask*lme_mask
            olme_mask= np.zeros_like(clme_mask)
            olme_mask1 = np.zeros_like(clme_mask)
            olme_mask2 = np.zeros_like(clme_mask)
            olme_mask3 = np.zeros_like(clme_mask)
            olme_mask4 = np.zeros_like(clme_mask)

            #points next to clme points that are deeper than Depth_lim or not in this lme
            olme_mask1[:,:-1,:]=clme_mask[:,1:,:] * np.logical_or(Dmasko[:,:-1,:],np.logical_not(clme_mask[:,:-1,:]))
            olme_mask2[:,1:,:] =clme_mask[:,:-1,:] * np.logical_or(Dmasko[:,1:,:],np.logical_not(clme_mask[:,1:,:]))
            olme_mask3[:,:,:-1]=clme_mask[:,:,1:] * np.logical_or(Dmasko[:,:,:-1],np.logical_not(clme_mask[:,:,:-1]))
            olme_mask4[:,:,1:] =clme_mask[:,:,:-1] * np.logical_or(Dmasko[:,:,1:],np.logical_not(clme_mask[:,:,1:]))
            olme_mask = olme_mask1 + olme_mask2 + olme_mask3 + olme_mask4
            olme_mask[olme_mask>1] = 1
            olme_mask=olme_mask*mask
            areao = np.sum(np.sum(olme_mask*DX*DY,axis=2),axis=1)





            #%%
            if not  'bgc' in EXPNAM:

                print('calculate pea')
                t0 = time.time()

                Zd_mask, _, _ = nemo.calculate_vertical_mask(200.)
                strat = coast.GriddedStratification(nemo)
                strat.calc_pea(nemo, Zd_mask)
                nemo.dataset['pea']=strat.dataset['PEA']
                t1 = time.time()
                print(t1 - t0)
            if 'bgc' in EXPNAM:
                bottom = nemo.dataset.variables['bottom_level'] -1
                #bottom[np.where(bottom<0.)] =0.
                #for j in range(ny):
                #    for i in range(xy):
                #        O2[:,j,i]=nemo.dataset['O2_o'][:,bottom[j,i],j,i]
                O2_bot = nemo.dataset['O2_o'].isel(z_dim = bottom)
                nemo.dataset['O2_bot'] = O2_bot
            #%%
            print('Processing indicators')
            if 'bgc' in EXPNAM:
                vars2=vars_bgc
                varname2=varname_bgc
            else:
                vars2=vars
                varname2=varname

            for ivar,var in enumerate(vars2):
                t0 = time.time()

                mask=clme_mask
                Area=area
                if 'o' in var:
                    mask=olme_mask
                    Area = areao
                print(varname2[ivar],vars2[ivar])

                if len(nemo.dataset[varname2[ivar]].sizes) == 4:
                    #data = nemo.dataset[varname2[ivar]].values[:,0,:,:]
                     data = mask*DX*DY*nemo.dataset[varname2[ivar]].isel(z_dim=0)
                else:
                     data = mask*DX*DY*nemo.dataset[varname2[ivar]]
                Var = (data.sum(dim='x_dim').sum(dim='y_dim') / Area)

                it=np.arange((year-ystart0)*12,(year-ystart0)*12+12)
                indicators[var][it] = Var.values

                print('data loaded')

                t1=time.time()
                print(t1-t0)

    #%%
    domain_outpath='/home/users/jholt/work/SENEMO/SENEMO_FC/'
    outname = f"{domain_outpath}physical_indicators_bgc_{LME_name}.p"
    with open(outname,'wb' ) as f:
        pickle.dump(indicators, f)
        f.flush()
#%%
#inname=f"{domain_outpath}physical_indicators.p"
#with open(inname,'rb' ) as f:

#   A['slmean']=slmean
#   A=pickle.load(f)
