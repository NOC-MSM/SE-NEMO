## RUNS

Summary of runs carried out as part of the SE-NEMO development

| RUN ID      | Coordinate | Mixing |  Tide | Rivers | TDISS | Wave Drag | HPG Scheme | Bathymetry | Length |
| :---        |    ----:   |   ---: | :---: |   ---: | :---: | :---:     | ---:       | :---:      | ---:   |
| Header      | Title      |        |       |        |       |           |            |            |        |



Variations in the runs are controlled in the header of the [`runscript.slurm`](https://github.com/NOC-MSM/SE-NEMO/blob/master/EXPREF/runscript.slurm).

**ZPS_REF_NOTIDE (a.k.a GS1p0_notide):**

```
#################### nemo runscript options ############################
# For info on the parameters see namelist_ref                          #
########################################################################
rn_rdt=600          ; ln_zps='.true.'      ; ln_tmx_itf='.true.'
ln_bt_auto='.true.' ; rn_bt_cmax=0.8        ; nn_baro=30
nn_mxlice=3         ; nn_z0_ice=1           ; ln_rnf_new='.false.'
ln_rstdate='.true.' ; ln_shlat2d='.false.'   ; nn_diaharm=1981
rn_Cd0=1.0e-3       ; ln_loglayer='.false.' ; ln_tide='.false.'
ln_boost='.true.'   ; ln_gls='.false.'       ; ln_int_wave_drag='.false.'
ln_hpg_djc='.false.' ;
########################################################################
```

**ZPS_REF_TIDE (a.k.a GS1p1_tide):**
```
#################### nemo runscript options ############################
# For info on the parameters see namelist_ref                          #
########################################################################
rn_rdt=600          ; ln_zps='.true.'      ; ln_tmx_itf='.false.'
ln_bt_auto='.true.' ; rn_bt_cmax=0.8        ; nn_baro=30
nn_mxlice=3         ; nn_z0_ice=1           ; ln_rnf_new='.false.'
ln_rstdate='.true.' ; ln_shlat2d='.false.'   ; nn_diaharm=1981
rn_Cd0=2.5e-3       ; ln_loglayer='.false.' ; ln_tide='.true.'
ln_boost='.true.'   ; ln_gls='.true.'       ; ln_int_wave_drag='.true.'
ln_hpg_djc='.false.' ;
########################################################################
```
N.B. GS1p1_tide uses the namelist_ref wave drag file, defined as:
```
cn_int_wave_drag = './INPUTS/tdiss_R025.nc'  ! filename for internal wave drag dissipation
```
and is not modified further in namelist_cfg_template.



**EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2 (a.k.a GS1p2_full):**

```
#################### nemo runscript options ############################
# For info on the parameters see namelist_ref                          #
########################################################################
rn_rdt=600          ; ln_zps='.false.'      ; ln_tmx_itf='.false.'
ln_bt_auto='.true.' ; rn_bt_cmax=0.8        ; nn_baro=30
nn_mxlice=3         ; nn_z0_ice=1           ; ln_rnf_new='.true.'
ln_rstdate='.true.' ; ln_shlat2d='.true.'   ; nn_diaharm=1981
rn_Cd0=2.5e-3       ; ln_loglayer='.true.'  ; ln_tide='.true.'
ln_boost='.false.'   ; ln_gls='.true.'      ; ln_int_wave_drag='.true.'
ln_hpg_djc='.true.' ;
########################################################################
```
N.B. In namelist_cfg_template, the wave drag file is defined for GS1p2_full as:
```
cn_int_wave_drag = 'INPUTS/tdiss_R025_fac2.nc',
```

Data can be found on JASMIN under `/gws/nopw/j04/class_vol2/senemo`.
