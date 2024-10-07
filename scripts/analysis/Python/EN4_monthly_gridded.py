import numpy as np
import matplotlib.pylab as plt
import xarray as xr
import coast
LME_Data=np.load('../Data/LME_gridinfo_V4.npz')
ystart = 1978
ystop = 2019
J_offset=186 #account for extra rows in eORCA if data is made for normal ORCA
fn_profile_config='../Config/example_en4_profiles.json'
in_path = '/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT/EN4.2.1/1978-2019'
out_path = '/projectsa/NEMO/jholt/SE-NEMO/ASSESSMENT/EN4.2.1/eORCA025_1978-2019'

fn_domain='/projectsa/NEMO/jholt/SE-NEMO/INPUTS/eORCA025_bathy_meter.nc'
for iLME in [57]:#range(66): # [31]:
    LME_Name= LME_Data["DOMNAM"][iLME]
    print(LME_Name)
    i_min = LME_Data['i_min'][iLME]
    i_max = LME_Data['i_max'][iLME]
    j_min = LME_Data['j_min'][iLME] + J_offset
    j_max = LME_Data['j_max'][iLME] + J_offset

    limits = [j_min, j_max, i_min, i_max]


    nemo = coast.Gridded(fn_data=fn_domain, config='example_nemo_grid_t.json')
    fname_in  = '{0}/{1}_{2}_{3}_EN4_PEA_SST_SSS_v1.nc'.format(in_path, LME_Name, ystart, ystop)
    fname_out = '{0}/{1}_{2}_{3}_EN4_eORCA025_PEA_SST_SSS_v1.nc'.format(out_path, LME_Name, ystart, ystop)
    isfile=True
    try:
        ds=xr.open_dataset(fname_in)
    except:
        isfile = False
    if isfile:
        pa=coast.Profile(dataset = ds)
        pa.grid_vars_mnth(nemo, ['sst', 'sss', 'pea'],limits=limits)
        nemo.dataset.to_netcdf(fname_out)