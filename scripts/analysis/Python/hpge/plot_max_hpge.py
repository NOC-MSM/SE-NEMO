#!/usr/bin/env python

import numpy as np
import xarray as xr
from utils import compute_masks
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as feature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy as ctp
from matplotlib import colors as c

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

# Parameters ------------------

sp1_vmin=0
sp1_vmax=2.

#set default font size
mpl.rcParams.update({'font.size': 16})

#set figure title font size
tfs=20

clonfixed=-80 #fixed central longitude for our preferred version of the Robinson projection
fixedprojection=ccrs.Robinson(central_longitude=clonfixed)

# Preferred colormaps
cmapsequential = truncate_colormap(plt.get_cmap('nipy_spectral'), minval=0.0, maxval=0.9, n=100)

#------------------------------------

hpge   = "/scratch/dbruciaf/SE-NEMO/hpge/u-cg602_hpge_se-nemo_r018-010-010_glo-r018-010_ant_opt_v3_3months_traldfoff/maximum_hpge_test.nc"
domcfg = "/data/users/dbruciaf/SE-NEMO/se-orca025/se-nemo-domain_cfg/domain_cfg_MEs.nc"

# Dealing with coordinates

ds_dom = xr.open_dataset(domcfg).squeeze()

# Computing land-sea masks
ds_dom = compute_masks(ds_dom, merge=True)
ds_dom = ds_dom.set_coords(["nav_lon", "nav_lat"])

nav_lon = ds_dom.nav_lon
nav_lat = ds_dom.nav_lat
tmask2D = ds_dom.tmask[0,:,:]

del ds_dom

ds_hpge = xr.open_dataset(hpge).squeeze()
ds_hpge = ds_hpge.assign_coords(nav_lon=nav_lon, nav_lat=nav_lat)
# Add x,y coords - fundamental for getting rid of discontinuity on lon grid
ds_hpge.coords["x"] = range(ds_hpge.dims["x"])
ds_hpge.coords["y"] = range(ds_hpge.dims["y"])
max_vel = ds_hpge.max_hpge_1
max_vel = max_vel.where(tmask2D == 1)

# Get rid of discontinuity on lon grid 
#(from https://gist.github.com/willirath/fbfd21f90d7f2a466a9e47041a0cee64)
# see also https://docs.xarray.dev/en/stable/examples/multidimensional-coords.html
after_discont = ~(max_vel.coords["nav_lon"].diff("x", label="upper") > 0).cumprod("x").astype(bool)
max_vel.coords["nav_lon"] = (
    max_vel.coords["nav_lon"]
    + 360 * after_discont
)

max_vel = max_vel.isel(x=slice(1, -1), y=slice(None, -1))

del ds_hpge

# Make the plot ---------------------------

def main():
    fig = plt.figure(figsize=(21, 15))
    #fig.suptitle('Max spurious currents', fontsize=tfs)
    plt.subplots_adjust(wspace=0.1,hspace=0.1)

    ax = fig.add_subplot(1, 1, 1, projection=fixedprojection)
    #ax.set_title(subplottitle1)
    ax.set_global()
    ax.set_aspect(1.2)

    #sets the land colour as grey (must be missing data)
    ax.set_facecolor("grey") 
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), 
                      draw_labels=True, 
                      xlocs=np.linspace(-160,200,19), 
                      ylocs=np.linspace(-80,80,17),
                      linewidth=2, 
                      color='gray', 
                      alpha=0.5, 
                      linestyle='--'
                     )
    gl.xlabels_top = False
    gl.ylabels_right = False

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Plot scalar field, here using pcolormesh.  Could substitute another plotting method, 
    # but please attempt to retain this choice of colormap

    
    #scalar_subplot = max_vel.plot.pcolormesh(
    #                 x="nav_lon", y="nav_lat",
    #                 ax=ax,
    #                 transform=ccrs.PlateCarree(), 
    #                 cmap=cmapsequential, 
    #                 vmin=sp1_vmin, 
    #                 vmax=sp1_vmax
    #)   
    pcol = ax.pcolormesh(
                        max_vel.nav_lon, 
                        max_vel.nav_lat,
                        max_vel*100.,
                        transform=ccrs.PlateCarree(), 
                        cmap=cmapsequential, 
                        vmin=sp1_vmin, 
                        vmax=sp1_vmax
    )   

    cbar_label = r"[$cm \, s^{-1}$]"
    cbaxes = inset_axes(ax, width="22%", height="4%", loc=6, borderpad=17)
    #cb = plt.colorbar(pcol, cax=cbaxes, ticks=[0., 2.5, 5.],orientation='horizontal', extend=cbar_extend)
    cb = plt.colorbar(pcol, cax=cbaxes, ticks=[0., 1., 2.],orientation='horizontal', extend='max')
    cb.ax.tick_params(color='white', labelcolor='white', labelsize=25, width=4, length=28, direction='in')
    cb.outline.set_edgecolor('white')
    cb.set_label(label=cbar_label,size=25, color='white')



    name = 'hpge_map.png'
    plt.savefig(name,bbox_inches="tight")


if __name__ == '__main__':
    main()

