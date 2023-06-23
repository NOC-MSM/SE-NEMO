# Shelf Enabled Global NEMO (SE-NEMO): MPI variations

The setup script has been tested and will checkout, compile and run the ORCA025 (NEMO 4.0.4) code on ARCHER2 for three versions of MPI (and two compilers): Cray-MPICH, GNU-MPICH, GNU-MPICH4 and GNU-openMPI.

**NB on this branch GNU-MPICH is the only tested (and currently working) version**

## Quick Start:

```
git clone git@github.com:NOC-MSM/SE-NEMO.git SE-NEMO-MPI-UPDATE

cd SE-NEMO-MPI-UPDATE
git checkout mpi_update

MPI_OPT='mpich' # other options are mpich4 | ompi

cd ../

./SE-NEMO-MPI-UPDATE/scripts/setup/se-eORCA025_setup -w $PWD/test_deploy_$MPI_OPT -x $PWD/test_deploy_$MPI_OPT -s $PWD/SE-NEMO-MPI-UPDATE -m archer2 -a $MPI_OPT

cd $PWD/test_deploy_$MPI_OPT/nemo/cfgs/se-eORCA025/
cp -rP EXPREF EXP_MYRUN
cd EXP_MYRUN
ln -s INPUTS/domain_cfg_mes.nc domain_cfg.nc # terrain following case (MES)
# alternativly use INPUTS/domain_cfg_zps.nc for the geopotential vertical coordinate (ZPS) case
```
Edit the project code and options in  `runscript.[slurm|mpirun]` then:
```
sbatch runscript.[slurm|mpirun] # the openMPI version can use either script
```
This will produce a 5 day mean output from the beginning of 1976. The run should take 15 minutes to complete once in the machine.

If using the ZPS case the following need to be set in the `runscript.[slurm|mpirun]`:
```
ln_zps='.true.'
ln_hpg_djc='.false.'
ln_loglayer='.false.'
ln_boost='.true.'
```

In the header of each runscript file there are three pre-defined core placements, for example in the `runscript.slurm` file:
```
#################### nemo runscript options ############################
# The following options are for the MES domain_cfg.nc                  #
########################################################################
if [ $SLURM_NNODES -eq 19 ]
then
   SRUN_CMD=`./slurm_setup -S 16 -s 16 -m 1 -C 1524 -g 4 -N 128 -t 01:00:00 -a n01-CLASS -j se-nemo`
elif [ $SLURM_NNODES -eq 29 ]
then
   SRUN_CMD=`./slurm_setup -S 16 -s 16 -m 1 -C 2973 -g 8 -N 128 -t 01:00:00 -a n01-CLASS -j se-nemo`
elif [ $SLURM_NNODES -eq 52 ]
then
   SRUN_CMD=`./slurm_setup -S 16 -s 16 -m 1 -C 6270 -g 9000 -N 128 -t 01:00:00 -a n01-CLASS -j se-nemo`
elif [ $SLURM_NNODES -eq 68 ]
then
   SRUN_CMD=`./slurm_setup -S 32 -s 16 -m 2 -C 8090 -g 9000 -N 128 -t 01:00:00 -a n01-CLASS -j se-nemo`
else
   exit
fi
########################################################################
```
**NB these pre-defined core placements are for the MES coordinates. Below this there are examples should you wish to run the ZPS coordinate configuration.**

So all that is required to switch between these is to alter `#SBATCH --nodes=XX` at the top of the script. By following this simple syntax, additional experiments can easily be added. To get information about about the use of `build_rankfile` just issue: `build_rankfile -h`. Likewise for the `runscript.slurm`: `slurm_setup -h`.

### Forcing data:

[SE-ORCA025](http://gws-access.ceda.ac.uk/public/jmmp_collab/)

_this is automatically transferred when the setup script is executed_

For ARCHER2 users these data are held under `/work/n01/shared/senemo` and `/work/n01/shared/nemo/FORCING` and are linked during the setup.

### Current progress:

At higher core counts it has become increasingly problematic to run NEMO/XIOS. The thought is that it may be linked to the default MPICH installation, so other MPI setups are being explored. The difficulties may also be attributed to balancing the NEMO to XIOS core count and memory available. Below is a summary of current progress.

_core count vs MPI_

|  MES    | GNU-MPICH                      | GNU-MPICH4                           | GNU-OMPI-srun|  GNU-OMPI-mpirun| Cray-MPICH |
| :----:  |  :----:                    |   :----:                         |:----:  |:----:  |:----:  |
| 1524    | 2.9 hrs/yr | | | |
| 2973    | 2.1 hrs/yr | | | |
| 6270    | 1.5 hrs/yr | | | |
| 8090    | 1.5 hrs/yr | | | |

|  ZPS    | GNU-MPICH                      | GNU-MPICH4                           | GNU-OMPI-srun|  GNU-OMPI-mpirun| Cray-MPICH |
| :----:  |  :----:                    |   :----:                         |:----:  |:----:  |:----:  |
| 1524    |  | | | |
| 2973    |  | | | |
| 6270    |  | | | |
| 8090    |  | | | |



### Note on iodef.xml

There are several *levers* you can pull with regards to XIOS management. One option is to increase the number of XIOS servers employed. Note that in NEMO version < 4.2 you must make sure that no XIOS server occupies land-only regions of the configuration. For example if I choose 8 XIOS servers then the SE-NEMO is divided up into 8 zonal `xios_server.exe` strips each ~152 points in *height* (i.e. 1221 j-points divided by 8). This decomposition works for the SE-NEMO configuration as the first XIOS server will have a small region of the Weddell Sea. If 16 XIOS servers are requested this becomes an issue - so single grid cell wet points are added to the domain_cfg file in Antarctica to overcome this). Other options when running into dificulties with XIOS can be adjusted in the `iodef.xml` file. The buffer size is optimised for performance by default. This can be either optimised for memory or can be scaled up or even both by adding the following code:

```
<variable_definition>
          <variable_group id="buffer">
             <variable id="optimal_buffer_size" type="string">memory</variable>
             <variable id="buffer_size_factor" type="double">4.0</variable>
          </variable_group>
          ...
```
          
