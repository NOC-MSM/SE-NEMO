## Generate bathymetric roughness field, following the basic approach of Jayne and St. Laurent (2001)

### Author: Chris Wilson, with input from Jason Holt and Jeff Polton.  9 Aug. 2023

**This is based on a previous estimate of David Byrne, but rewritten from scratch, with minor changes.**

As it’s quite a large computation and because the JupyterLab (JASMIN Notebook Service) seems slightly unstable for long jobs, there is both a Jupyter notebook and a plain python script here.  The latter was used on a single node, multi-core workstation (the JASMIN Lotus cluster) to do the calculation.

The algorithm involves the following:

* Read the bathymetry data (assumed to be on a grid with much finer gridscale than the NEMO model.
* Read the NEMO model grid (via its domain_cfg file).
* Use the Arakawa C-grid layout of NEMO and that its bathymetry is defined on T-points to define the lon-lat bounds of each of its gridcells using the relevant u- and v-point lons and lats in that cell and adjacent cells.
* Make a mask of the NEMO ocean(1)/land(0) T-grid cells, based on the “top level” variable.
* Loop over each of the NEMO ocean cells to:
  * find all the bathymetry data within the lon-lat boundaries of the cell (this is an approximation, as the NEMO grid (e.g. tripolar ORCA) may not necessarily align perfectly with the lon-lat grid.
  * Using nonlinear least-squares optimisation, fit the surface H=a+bx+cy+dxy to the bathymetry, by solving for a, b, c, d.   Note that the coordinates x,y here relate to the ii, jj fine-grid indices of the bathymetry grid “patch”.  As the patch is assumed to be very small, these ii, jj approximate to Great Circle distances x,y anyway.    The surface fitted is a hyperbolic paraboloid.  See example in Jupyter notebook.  For small d, this approximates a bilinear fit of a plane to the bathymetry of each patch.
  * Define the anomaly from the fitted surface, h’=bathymetry-H,  and calculate the root-mean-square h’ defined over all the ii,jj of the patch, as the **roughness, h**.

This python script can be run to calculate the bathymetric roughness:
```
barebones_calc_bathymetric_roughness.py
```
and this Jupyter notebook is similar to the above, but includes test plots and narrative:
```
calc_bathymetric_roughness.ipynb
```

Notes:

1. It’s possible to use different choices of bathymetry data and model grid - just edit the script.  Here, we used GEBCO2023, which is 15’ lat-lon resolution, and the SENEMO grid corresponding to the domain_cfg for GS1p2_full, which is eORCA025 global, tripolar.
2. It might be possible to speed up the calculation further, using Dask, but I struggled to get it to work.  The need to iterate over the NEMO model grid, plus to optimally fit the function H and calculate h seemed to cause headaches.  

Acknowledgement:
**GEBCO Compilation Group (2023) GEBCO 2023 Grid (doi:10.5285/f98b053b-0cbc-6c23-e053-6c86abc0af7b)**
