#!/bin/bash
#SBATCH --job-name=SE-ZIP
#SBATCH --time=24:00:00
#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

. /work/n01/n01/jdha/modules/gnu-mpich

if [ ! -d OUTPUTS_ZIP ]
then
	mkdir OUTPUTS_ZIP
fi

if [ ! -d OUTPUTS_PROCESSED ]
then
	mkdir OUTPUTS_PROCESSED
fi

cd OUTPUTS

for i in *
do
	nccopy -4 -d 1 -c x/64,y/64 $i ../OUTPUTS_ZIP/$i
	mv $i ../OUTPUTS_PROCESSED/
done

