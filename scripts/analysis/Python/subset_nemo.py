#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 17:05:04 2023

@author: jholt
"""

coast_dir='/home/users/jholt/Git/COAsT/'
import sys
sys.path.insert(0,coast_dir)
import matplotlib.pylab as plt
import coast
import numpy as np
import xarray as xr
import scipy.io
import os
os.system("module load jasmin-sci")
names,dpaths,DOMS,_  = coast.experiments(experiments='experiments.json')

nemo_t=coast.Gridded(fn_domain=DOMS[0],config='example_nemo_grid_t.json',no_depths=True)


x_min=-60;x_max=-40;y_min=-60;y_max=-40

j,i,_=nemo_t.find_j_i_list(lon=[x_min,x_max,x_max,x_min],lat=[y_min,y_min,y_max,y_max])
imin=min(i)
imax=max(i)
jmin=min(j)
jmax=max(j)

nemo_t.subset(y_dim=range(jmin,jmax),x_dim=range(imin,imax))

ystart=2010
ystop=2019

fn_nemo_dat_t= coast.nemo_filename_maker(dpaths[0],ystart,ystop,grid='T')

outdir=f"/home/users/jholt/SENEMO/Subsets_{imin}_{imax}_{jmin}_{jmax}"
if not os.path.exists(outdir):
    os.mkdir(outdir)
for name in fn_nemo_dat_t:
    name_out= f"{outdir}/{os.path.basename(name)}"
    os.system(
        f"ncea -O -d x,{imin},{imax} -d y,{jmin},{jmax} {name} {name_out}"
        )