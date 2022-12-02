## RUNS

Summary of runs carried out as part of the SE-NEMO development

| RUN ID      | Coordinate | Mixing |  Tide | Rivers | TDISS | Wave Drag | HPG Scheme | Bathymetry | Length |
| :---        |    ----:   |   ---: | :---: |   ---: | :---: | :---:     | ---:       | :---:      | ---:   |
| Header      | Title      |        |       |        |       |           |            |            |        |



Variations in the runs are controlled in the header of the [`runscript.slurm`](https://github.com/NOC-MSM/SE-NEMO/blob/master/EXPREF/runscript.slurm).

```
#################### nemo runscript options ############################
# For info on the parameters see namelist_ref                          #
########################################################################
rn_rdt=600          ; ln_zps='.false.'      ; ln_tmx_itf='.false.'
ln_bt_auto='.true.' ; rn_bt_cmax=0.8        ; nn_baro=30
nn_mxlice=3         ; nn_z0_ice=1           ; ln_rnf_new='.false.'
ln_rstdate='.true.' ; ln_shlat2d='.true.'   ; nn_diaharm=1981
rn_Cd0=2.5e-3       ; ln_loglayer='.false.' ; ln_tide='.true.'
ln_boost='.true.'   ; ln_gls='.true.'       ; ln_int_wave_drag='.true.' 
ln_hpg_djc='.true.' ; cn_int_wave_drag = './INPUTS/tdiss_R025.nc'
########################################################################
```

Data can be found on JASMIN under `/gws/nopw/j04/class_vol2/senemo`.
