#!/bin/bash
#SBATCH --job-name=oCMPR
#SBATCH --output=out-CMPR-OUT-%j.out
#SBATCH --time=06:00:00
#SBATCH --account=n01-pml
#SBATCH --partition=serial
#SBATCH --qos=serial

# sbatch --export=ARCHIVEDIR=$ARCHIVEDIR,y0=$y0,m0=$m0 compress_output.slurm
# nccopy options: -s: shuffling, -u: removes unlimited dimension

module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel

INDIR=$PWD

mkdir -p $ARCHIVEDIR

rm timing_compr_out.txt
touch timing_compr_out.txt
echo start: >> timing_compr_out.txt
ls -lhrt timing_compr_out.txt >> timing_compr_out.txt

for file in SENEMO_S?_1?_$y0$m0str*.nc 
#for file in SENEMO_1?_$y0$m0str*.nc 
do
  echo compressing $file
  nccopy -s -u -k 4 -d 5 $INDIR/$file $ARCHIVEDIR/$file && rm $INDIR/$file || echo $file compression failed !!!!
done

touch timing_compr_out.txt
echo end: >> timing_compr_out.txt
ls -lhrt timing_compr_out.txt >> timing_compr_out.txt
