# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import xarray as xr
import matplotlib.pylab as plt
import numpy as np

fn1="/work/n01/n01/cwi/projects/senemo/buildIWD160823/testIWD/nemo/cfgs/se-eORCA025/EXPIWD02/OUTPUTS/SENEMO_1y_19780101_19781231_grid_T_2D.nc"
f1=xr.open_dataset(fn1)
#fn2="/work/n01/n01/jdha/temp_se_nemo/nemo/cfgs/se-eORCA025/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/OUTPUTS_PROCESSED/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
#f2=xr.open_dataset(fn2)

fn3="/work/n01/n01/cwi/projects/senemo/buildIWD160823/testIWD/nemo/cfgs/se-eORCA025/EXPIWD01/OUTPUTS/SENEMO_1y_19780101_19781231_grid_T_2D.nc"
f3=xr.open_dataset(fn3)
fn4="/home/n01/n01/jholt/work/SENEMO//nemo/cfgs/se-eORCA025/EXP_GS1p6_full_IWD_soenhance//OUTPUTS/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
f4=xr.open_dataset(fn4)


#fn5="/work/n01/n01/cwi/projects/senemo/buildIWD160823/testIWD/nemo/cfgs/se-eORCA025/EXPIWD05/OUTPUTS/SENEMO_1y_19780101_19781231_grid_T_2D.nc"
#f5=xr.open_dataset(fn5)

fn5="/work/n01/n01/jdha/scratch/SN_JRA_RIV/nemo/cfgs/se-eORCA025/EXP_TEST/OUTPUTS_ZIP/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
f5=xr.open_dataset(fn5)




M2_1=np.sqrt(f1.M2x.values[:,:]**2+f1.M2y.values[:,:]**2)

#M2_2=np.sqrt(f2.M2x.values[:,:]**2+f2.M2y.values[:,:]**2)

M2_3=np.sqrt(f3.M2x.values[:,:]**2+f3.M2y.values[:,:]**2)

M2_4=np.sqrt(f4.M2x.values[:,:]**2+f4.M2y.values[:,:]**2)

M2_5=np.sqrt(f5.M2x.values[:,:]**2+f5.M2y.values[:,:]**2)