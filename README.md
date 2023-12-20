# Shelf Enabled Global NEMO (SE-NEMO)

The setup script has been tested and will checkout, compile and run the ORCA025 (NEMO 4.0.4 or 4.2.1) code on: ARCHER2 for Cray-MPICH and GNU-MPICH, and Anemone for iFort. 

Configuration files for SE-NEMO project

Base Configuration: GO8p6 at NEMO 4.0.4 or NEMO 4.2.1

## Quick Start:
On ARCHER2
```
git clone git@github.com:NOC-MSM/SE-NEMO.git
./SE-NEMO/scripts/setup/se-eORCA025_setup -p $PWD/test  -r $PWD/cfgs/SE-NEMO -n 4.0.4 -x 2.5 -m archer2 -a mpich -c gnu
# or for 4.2.1
#./SE-NEMO/scripts/setup/se-eORCA025_setup -p $PWD/test  -r $PWD/cfgs/SE-NEMO -n 4.2.1 -x 2 -m archer2 -a mpich -c gnu
cd test/nemo/cfgs/se-eORCA025/
cp -rP EXPREF EXP_MYRUN
cd EXP_MYRUN
ln -s ../INPUTS/domain_cfg_zps.nc domain_cfg.nc # ZPS or domain_cfg_mes_v2.nc MES
# or for 4.2.1
# ln -s ../INPUTS/domain_cfg_zps_nohalos.nc domain_cfg.nc # ZPS or domain_cfg_mes_nohalos.nc MES

```
or if using ANEMONE, replace use options:
```
-m anemone -a impi -c ifort
```
Edit the project code and options in  `runscript.slurm` then:
```
sbatch runscript.slurm
```
This will produce a 5 day mean output from the beginning of 1976. The run should take 15 minutes to complete once in the machine.

### Forcing data:

[SE-ORCA025](http://gws-access.ceda.ac.uk/public/jmmp_collab/)

_this is automatically transferred when the setup script is executed_

For ARCHER2 users these data are held under `/work/n01/shared/senemo` and `/work/n01/shared/nemo/FORCING` and are linked during the setup.

### Runs:

A list of ongoing simluations can be found [here](https://github.com/NOC-MSM/SE-NEMO/blob/master/RUNS.md)

### Outputs:

On JASMIN: /gws/nopw/j04/class_vol2/senemo

### Important:

For 4.0.4 code: the `MY_SRC_GO8_FROZEN` directory should only contain the re-ROSEd code from G08. If editing any of these routines, copy the file from `MY_SRC_GO8_FROZEN` to `MY_SRC` and edit it there. The build process will copy the contents of `MY_SRC_GO8_FROZEN` to `cfgs/se-nemo/MY_SRC` before copying `MY_SRC` to `cfgs/se-nemo/MY_SRC`. Any files with the same name will be overwritten.
