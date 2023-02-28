#!/bin/bash
#SBATCH --job-name=RBLD
#SBATCH --output=out-rbld-%j.out
#SBATCH --time=06:00:00
#SBATCH --account=n01-pml
##SBATCH --partition=serial
##SBATCH --qos=serial
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=2
##SBATCH --cpus-per-task=64
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
##SBATCH --mem=125G

# submit with:
# sbatch --export=ITER=4464,SPIN=1 rebuild_restart.slurm

#rsync -aP /work/n01/n01/gig/eORCA025/rebuild_nemo.exe .

#module load PrgEnv-gnu
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7

CURR_DIR=$PWD
SERIAL=$(printf "%08d" $ITER)
INDIR=$PWD/RESTARTS
OUTDIR=$PWD/RESTARTS_PROCESSED
ndom=$(ls -1 $INDIR/SENEMO_S"$SPIN"_"$SERIAL"_restart_trc_????.nc | wc -l)
#ndom=$(ls -1 $INDIR/SENEMO_"$SERIAL"_restart_trc_????.nc | wc -l)

mkdir -p $OUTDIR

rm nam_rebuild_trc
rm nam_rebuild_ice
rm nam_rebuild
#echo -e "&nam_rebuild\nfilebase='"$INDIR"/SENEMO_"$SERIAL"_restart_trc'\nndomain=$ndom\n/" >> nam_rebuild_trc
#echo -e "&nam_rebuild\nfilebase='"$INDIR"/SENEMO_"$SERIAL"_restart_ice'\nndomain=$ndom\n/" >> nam_rebuild_ice
#echo -e "&nam_rebuild\nfilebase='"$INDIR"/SENEMO_"$SERIAL"_restart'\nndomain=$ndom\n/" >> nam_rebuild

echo -e "&nam_rebuild\nfilebase='"$INDIR"/SENEMO_S"$SPIN"_"$SERIAL"_restart_trc'\nndomain=$ndom\n/" >> nam_rebuild_trc
echo -e "&nam_rebuild\nfilebase='"$INDIR"/SENEMO_S"$SPIN"_"$SERIAL"_restart_ice'\nndomain=$ndom\n/" >> nam_rebuild_ice
echo -e "&nam_rebuild\nfilebase='"$INDIR"/SENEMO_S"$SPIN"_"$SERIAL"_restart'\nndomain=$ndom\n/" >> nam_rebuild

echo $SERIAL

rm timing_rbld
touch timing_rbld
echo start: >> timing_rbld
ls -lhrt timing_rbld >> timing_rbld


# THIS RUNS!
# restaret_trc takes 5h 26m on compressed files, 2h on uncompressed files
echo ...REBUILDING RESTARTS...
srun --distribution=block:block --hint=nomultithread --unbuffered ./rebuild_nemo.exe nam_rebuild_trc
srun --distribution=block:block --hint=nomultithread --unbuffered ./rebuild_nemo.exe nam_rebuild_ice
srun --distribution=block:block --hint=nomultithread --unbuffered ./rebuild_nemo.exe nam_rebuild

touch timing_rbld
echo rebuild completed >> timing_rbld
ls -lhrt timing_rbld >> timing_rbld

echo Done.
