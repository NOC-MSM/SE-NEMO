#!/usr/bin/env python

#     |------------------------------------------------------------|
#     | This module creates a 2D field of maximum spurious current |
#     | in the vertical and in time after an HPGE test.            |
#     | The resulting file can be used then to optimise the rmax   |
#     | of Multi-Envelope vertical grids.                          |
#     |                                                            |
#     | Author: Diego Bruciaferri                                  |
#     | Date and place: 07-09-2021, Met Office, UK                 |
#     |------------------------------------------------------------|


import os
from os.path import join, isfile, basename, splitext
import glob
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import xarray as xr
from utils import compute_masks

# ==============================================================================
# Input files
# ==============================================================================

# Folder path containing HPGE spurious currents velocity files 
MAINdir = '/scratch/dbruciaf/SE-NEMO/hpge/'
HPGElst = 'u-cg602_hpge_se-nemo_r018-010-010_glo-r018-010_ant_opt_v3_3months_traldfoff'
FileVel = 'hpge_timeseries.nc'

# ==============================================================================
# OPENING fig
fig, ax = plt.subplots(figsize=(16,9))

# Loading timeseries
ds  = xr.open_dataset(join(join(MAINdir,HPGElst),FileVel)).squeeze()

ax.plot(np.arange(1,32), ds.max_u*100., linestyle="-", linewidth=5, color='blue', label='max{$| \mathbf{u} |$}')
ax.plot(np.arange(1,32), ds.u_99p_loc*100., linestyle="--", linewidth=5, color='gold', label='99%{$| \mathbf{u} |$}')
ax.plot(np.arange(1,32), ds.avg_u*100., linestyle="--", linewidth=5, color='red', label='$V_L^{-1} \int_{V_L} | \mathbf{u} | \mathrm{d}V$')
       
plt.rc('legend', **{'fontsize':35})
ax.legend(loc=5, ncol=1, frameon=False)
ax.set_xlabel('Days', fontsize=35)
ax.set_ylabel('[$cm\;s^{-1}$]', fontsize=35)
ax.tick_params(axis='both',which='major', labelsize=30)
ax.set_xlim(1.,32)
ax.set_ylim(0.,6.)
ax.grid(True)
name = 'hpge_timeseries.png'
plt.savefig(name, bbox_inches="tight")
