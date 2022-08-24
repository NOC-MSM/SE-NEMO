#!/bin/bash
# Command line arguments:
#   $1: rundir
#   $2: simulation year

#INPUT_PATH=/work/n01/n01/...
INDIR_Ndep=/work/n01/shared/yuti/se-ORCA025_Ndep/
INDIR_Fedep=/work/n01/shared/yuti/se-ORCA025_Fedep/
INDIR_pCO2a=/work/n01/shared/yuti/se-ORCA025_pCO2a/
INDIR_ADY=/work/n01/shared/gig/eORCA025_INPUTS/eORCA025_ADY_clim/

INDIR_JDHA=/work/n01/n01/jdha/scratch/eORCA1/nemo/cfgs/eORCA1/EXP01
INDIR_JDHA_JRA=/work/n01/n01/jdha/scratch/eORCA1/nemo/cfgs/eORCA1/EXP01/INPUTS/JRA/JRA_v1.5.0_rechunk

RUNPATH=$PWD
yn=$1             # curr year
yb=$(( $yn-1 ))   # previous year
ya=$(( $yn+1 ))   # next year

ln -sf  $INDIR_Ndep/eORCA025_N_dep_ISIMIP_y${yn}.nc $RUNPATH
ln -sf  $INDIR_Ndep/eORCA025_N_dep_ISIMIP_y${yb}.nc $RUNPATH
ln -sf  $INDIR_Ndep/eORCA025_N_dep_ISIMIP_y${ya}.nc $RUNPATH

ln -sf  $INDIR_Fedep/eORCA025_Fe_dep_GESAMP.nc $RUNPATH

ln -sf  $INDIR_pCO2a/eORCA025_v2_pCO2a_global_y${yn}.nc $RUNPATH
ln -sf  $INDIR_pCO2a/eORCA025_v2_pCO2a_global_y${yb}.nc $RUNPATH
ln -sf  $INDIR_pCO2a/eORCA025_v2_pCO2a_global_y${ya}.nc $RUNPATH

ln -sf  $INDIR_ADY/eORCA025-CCI-ady-01-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m01.nc
