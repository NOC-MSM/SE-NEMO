## RUNS

| RUN ID      | Coordinate | Mixing |  Tide | Rivers | TDISS | Wave Drag | HPG Scheme | Bathymetry | 
| :---        |    ----:   |   ---: | :---: |   ---: | :---: | :---:     | ---:       | :---:      |
| Header      | Title      |        |       |        |       |           |            |            |



```
#################### nemo runscript options ############################
# For info on the parameters see namelist_ref                          #
########################################################################
rn_rdt=600          ; ln_zps='.false.'      ;
ln_bt_auto='.true.' ; rn_bt_cmax=0.8        ; nn_baro=30
nn_mxlice=3         ; nn_z0_ice=1           ; ln_rnf_new='.false.'
ln_rstdate='.true.' ; ln_shlat2d='.true.'   ; nn_diaharm=1981
rn_Cd0=2.5e-3       ; ln_loglayer='.false.' ; ln_tide='.true.'
ln_boost='.true.'   ; ln_gls='.true.'       ; ln_int_wave_drag='.false.'
########################################################################
```
