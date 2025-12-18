

import numpy as np
import scipy.io
import matplotlib.pylab as plt
import time
import xarray as xr
import pickle
import sys


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
vars = ['mld','sst','sss','ssso','nitrate']
varname =['mldr10_1','temperature','salinity','salinity','N3_n']

vars_bgc = ['nitrate']
varname_bgc =['N3_n']

#%%

ilme=21
nlme=66
ntmax=(max(year_stop)-min(year_start)+1)*12

#%%
for ilme in [6]:#,23]:
    indicators = xr.DataArray(dims='t_dim')
    for var in vars:
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
    for i in [0,2]:#1,2,3]:#, 1, 2, 3]:
        EXPNAM = names[i]
        ystart = year_start[i]
        ystop = year_stop[i]
        for year in range(year_start[i],year_stop[i]):
            print('opening',EXPNAM)
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


            for year in range(ystart, ystop + 1):
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
            nemos[i] = coast.Gridded(fn_data=fn_nemo_dat, fn_domain=fn_nemo_dom, config=fn_config_t_grid,
                                 multiple=True,no_depths=True,lims=lims)
            t1 = time.time()
            print(t1-t0)
        #%%
        print('Make blank object')
        nemo = coast.Gridded(config=fn_config_t_grid,no_depths=True)
        nemo_bgc = coast.Gridded(config=fn_config_t_grid, no_depths=True)

        #%%
        print('Concatenating')
        if 1 in list(nemos.keys()):
            nemo.dataset = xr.concat([nemos[0].dataset,nemos[1].dataset],dim='t_dim')
        else:
            nemo.dataset = nemos[0].dataset
        if 3 in list(nemos.keys()):
            nemo_bgc.dataset = xr.concat([nemos[2].dataset, nemos[3].dataset], dim='t_dim')
        else:
            nemo_bgc.dataset = nemos[2].dataset
        nemo.dataset = xr.merge([nemo.dataset,nemo_bgc.dataset])
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
        print('calculate pea')

        nyear=int(nemo.dataset.sizes["t_dim"]/12)
        t0 = time.time()
        Zd_mask, _,_ = nemo.calculate_vertical_mask(200.)
        strat = coast.GriddedStratification(nemo)
        first=True
        for iy in range(nyear):
            print("Calc pea", iy)
            it = np.arange((iy) * 12, (iy) * 12 + 12).astype(int)
            nemo2=nemo.subset_as_copy(t_dim=it)
            strat2 = coast.GriddedStratification(nemo2)
            strat2.calc_pea(nemo2, Zd_mask)
            if first :
                strat.dataset = strat2.dataset
                first=False
            else:
                strat.dataset=xr.concat([strat.dataset, strat2.dataset], dim='t_dim')
        nemo.dataset['pea']=strat.dataset['PEA']
        t1 = time.time()
        print(t1 - t0)

        #%%
        print('Processing indicators')

        indicators2 = xr.DataArray(dims='t_dim')

        for ivar,var in enumerate(vars):
            t0 = time.time()

            mask=clme_mask
            Area=area

            if 'o' in var:
                mask=olme_mask
                Area = areao
            print(varname[ivar],vars[ivar])

            if len(nemo.dataset[varname[ivar]].sizes) == 4:
                #data = nemo.dataset[varname[ivar]].values[:,0,:,:]
                 data = mask*DX*DY*nemo.dataset[varname[ivar]].isel(z_dim=0)
            else:
                 data = mask*DX*DY*nemo.dataset[varname[ivar]]
            Var = (data.sum(dim='x_dim').sum(dim='y_dim') / Area)
            indicators2[var][:ntimes]=Var.values
        indicators = xr.concat(indicators,indicators2, dim='t_dim')
            print('data loaded')

            t1=time.time()
            print(t1-t0)

    #%%
domain_outpath='/home/users/jholt/work/SENEMO/SENEMO_FC/'
outname = f"{domain_outpath}physical_indicators_{LME_name}.p"
with open(outname,'wb' ) as f:
    pickle.dump(indicators, f)
#%%
#inname=f"{domain_outpath}physical_indicators.p"
#with open(inname,'rb' ) as f:

#   A['slmean']=slmean
#   A=pickle.load(f)