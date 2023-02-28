#!/bin/bash

# bash clean_restarts.sh $m0
# deletes restarts from 2 months before
# must be adapted to count also for the year

#m0=$1
ITER=$1
INDIR=$PWD

echo cleaning restarts from previous two months, ITER: $ITER

#mm=$(( m0-2 ))
#if [ $mm2 -le 0 ]
#then
#  mm2=1
#fi
#mm2str=$(printf %02d $mm2)

#source nn_itend_${mm2str} #reads ITERmm2 value
SERIAL=$(printf "%08d" $ITER)

rm $INDIR/RESTARTS/SENEMO_S?_${SERIAL}_restart_????.nc
rm $INDIR/RESTARTS/SENEMO_S?_${SERIAL}_restart_ice_????.nc
rm $INDIR/RESTARTS/SENEMO_S?_${SERIAL}_restart_trc_????.nc

rm $INDIR/RESTARTS/SENEMO_${SERIAL}_restart_????.nc
rm $INDIR/RESTARTS/SENEMO_${SERIAL}_restart_ice_????.nc
rm $INDIR/RESTARTS/SENEMO_${SERIAL}_restart_trc_????.nc

