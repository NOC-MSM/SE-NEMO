# MY_SRC

<<<<<<< HEAD
Contains developement code for SE-NEMO
Square brackets show comment markers in the code to find and identify changes.

dynspg_ts.F90 - [davbyr] Applies internal wave drag prmtrsation to barotropic component.
=======
Contains development code for SE-NEMO
>>>>>>> master

dtatsd.F90 - vertical interpolation of ICs\
istate.F90 - veritcal interpolation of ICs\
par_oce.F90 - veritcal interpolation of ICs\
sbcrnf.F90 - bugfixes for missing input directory

sbctide.F90
   - [NB] Applies Love Number from namelist.
   - [NB] Applies long period tide forcing (previously set to zero).

tide.h90
   - [NB] Updated tidal potential forcing data. More constituents, Schureman method.

tide_mod.F90
   - [NB] Updated nodal factor equation.

tideini.F90
   - [NB] Reads Love number from namelist and outputs to Ocean.output.
   - [davbyr] Reads/stores ln_int_wave_drag (switches on internal wave drag for tides)
