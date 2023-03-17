#!/bin/bash

# bash clean_restarts.sh $m0 $nn_spin_cycle
# deletes restarts from 2 months before
# must be adapted to count also for the year

ITER=$1
SPIN=$2
INDIR=$PWD

echo cleaning restarts from previous two months, ITER: $ITER

SERIAL=$(printf "%08d" $ITER)

if [ $SPIN -ge 4 ] #NB, hardcoded 3y spin cycle
then
  rm $INDIR/RESTARTS/SENEMO_S?_${SERIAL}_restart_????.nc
  rm $INDIR/RESTARTS/SENEMO_S?_${SERIAL}_restart_ice_????.nc
  rm $INDIR/RESTARTS/SENEMO_S?_${SERIAL}_restart_trc_????.nc
else
  rm $INDIR/RESTARTS/SENEMO_${SERIAL}_restart_????.nc
  rm $INDIR/RESTARTS/SENEMO_${SERIAL}_restart_ice_????.nc
  rm $INDIR/RESTARTS/SENEMO_${SERIAL}_restart_trc_????.nc
fi
