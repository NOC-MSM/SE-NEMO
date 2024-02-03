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
ln_hpg_djc='.false.';                         ln_calc_tdiss='.false.'
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
ln_hpg_djc='.false.';                         ln_calc_tdiss='.false.'
########################################################################
```
N.B. GS1p1_tide does not calculate tdiss (internal wave drag) online, 
but also does not use the default input file 
in namelist_ref (cn_int_wave_drag = './INPUTS/tdiss_R025.nc').
Instead, that is modified by namelist_cfg_template to a field that's 
double in amplitude:
```
cn_int_wave_drag = 'INPUTS/tdiss_R025_fac2.nc',
```


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
ln_hpg_djc='.true.';                         ln_calc_tdiss='.true.'
########################################################################
```
N.B. GS1p2_full does calculate tdiss (internal wave drag) online.
It's necessary for both ln_int_wave_drag and ln_calc_tdiss to be '.true.'
to reach the appropriate subroutine.

Therefore, although there is a definition in namelist_cfg_template of 
```
cn_int_wave_drag = 'INPUTS/tdiss_R025_fac2.nc',
```
it is not used in this case.


Data can be found on JASMIN under `/gws/nopw/j04/class_vol2/senemo`.
