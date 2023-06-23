#!/bin/bash
# Command line arguments:
#   $1: rundir
#   $2: simulation year

#INPUT_PATH=/work/n01/n01/...
INDIR_Ndep=/work/n01/shared/yuti/se-ORCA025_Ndep/
INDIR_Fedep=/work/n01/shared/yuti/se-ORCA025_Fedep/
INDIR_pCO2a=/work/n01/shared/yuti/se-ORCA025_pCO2a/
INDIR_ADY=/work/n01/shared/yuti/se-ORCA025_ADY_clim/
INDIR_RIVERS=/work/n01/shared/yuti/se-ORCA025_JRA_BGC/dep_spread/

INDIR_JDHA=/work/n01/n01/jdha/scratch/eORCA1/nemo/cfgs/eORCA1/EXP01
INDIR_JDHA_JRA=/work/n01/n01/jdha/scratch/eORCA1/nemo/cfgs/eORCA1/EXP01/INPUTS/JRA/JRA_v1.5.0_rechunk

RUNPATH=$PWD
yn=$1             # curr year
yb=$(( $yn-1 ))   # previous year
ya=$(( $yn+1 ))   # next year

ln -sf  $INDIR_Ndep/eORCA025_r015-r010_007_004v2_N_dep_ISIMIP_y${yn}.nc $RUNPATH/eORCA025_N_dep_ISIMIP_y${yn}.nc
ln -sf  $INDIR_Ndep/eORCA025_r015-r010_007_004v2_N_dep_ISIMIP_y${yb}.nc $RUNPATH/eORCA025_N_dep_ISIMIP_y${yb}.nc
ln -sf  $INDIR_Ndep/eORCA025_r015-r010_007_004v2_N_dep_ISIMIP_y${ya}.nc $RUNPATH/eORCA025_N_dep_ISIMIP_y${ya}.nc

#cheat to do 2019 & 2020
#ln -sf  $INDIR_Ndep/eORCA025_r015-r010_007_004v2_N_dep_ISIMIP_y2018.nc $RUNPATH/eORCA025_N_dep_ISIMIP_y${yn}.nc
#ln -sf  $INDIR_Ndep/eORCA025_r015-r010_007_004v2_N_dep_ISIMIP_y2018.nc $RUNPATH/eORCA025_N_dep_ISIMIP_y${yb}.nc
#ln -sf  $INDIR_Ndep/eORCA025_r015-r010_007_004v2_N_dep_ISIMIP_y2018.nc $RUNPATH/eORCA025_N_dep_ISIMIP_y${ya}.nc

ln -sf  $INDIR_Fedep/eORCA025_r015-r010_007_004v2_Fe_dep_GESAMP.nc $RUNPATH/eORCA025_Fe_dep_GESAMP.nc

#ln -sf  $INDIR_RIVERS/ORCA025_rivers_Antar_Green_BGC_y${yn}.nc ${RUNPATH}/JRA_BGC_y${yn}.nc
#ln -sf  $INDIR_RIVERS/ORCA025_rivers_Antar_Green_BGC_y${ya}.nc ${RUNPATH}/JRA_BGC_y${ya}.nc
#ln -sf  $INDIR_RIVERS/ORCA025_rivers_Antar_Green_BGC_y${yb}.nc ${RUNPATH}/JRA_BGC_y${yb}.nc
ln -sf $INDIR_RIVERS/ORCA025_rivers_AG_BGC_dep_spread_y${yn}.nc ${RUNPATH}/JRA_BGC_y${yn}.nc
ln -sf $INDIR_RIVERS/ORCA025_rivers_AG_BGC_dep_spread_y${ya}.nc ${RUNPATH}/JRA_BGC_y${ya}.nc
ln -sf $INDIR_RIVERS/ORCA025_rivers_AG_BGC_dep_spread_y${yb}.nc ${RUNPATH}/JRA_BGC_y${yb}.nc
ln -sf  /work/n01/shared/se-eORCA025/eORCA_R025_runoff_v1.0.nc  ./runoff_1m_nomask.nc


ln -sf  $INDIR_pCO2a/eORCA025_r015-r010_007_004v2_pCO2a_global_y${yn}.nc $RUNPATH/eORCA025_pCO2a_global_y${yn}.nc
ln -sf  $INDIR_pCO2a/eORCA025_r015-r010_007_004v2_pCO2a_global_y${yb}.nc $RUNPATH/eORCA025_pCO2a_global_y${yb}.nc
ln -sf  $INDIR_pCO2a/eORCA025_r015-r010_007_004v2_pCO2a_global_y${ya}.nc $RUNPATH/eORCA025_pCO2a_global_y${ya}.nc

ln -sf  $INDIR_ADY/eORCA025_r015-r010_007_004v2_ady.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-01-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m01.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-02-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m02.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-03-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m03.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-04-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m04.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-05-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m05.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-06-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m06.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-07-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m07.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-08-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m08.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-09-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m09.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-10-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m10.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-11-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m11.nc
#ln -sf  $INDIR_ADY/eORCA025-CCI-ady-12-broadband-climatology.nc ${RUNPATH}/eORCA025-CCI-ady-broadband-climatology_m12.nc
