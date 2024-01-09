#!/bin/bash
#SBATCH --job-name=SENEMO
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=24


module load NEMO/prg-env

echo $PWD

export OMP_NUM_THREADS=1
export I_MPI_SHM=off
unset I_MPI_PMI_LIBRARY


#################### nemo runscript options ############################
# For info on the parameters see namelist_ref                          #
########################################################################
rn_rdt=600          ; ln_zps='.false.'      ; ln_tmx_itf='.false.'
ln_bt_auto='.true.' ; rn_bt_cmax=0.8        ; nn_baro=30
nn_mxlice=3         ; nn_z0_ice=1           ; ln_rnf_new='.true.'
ln_rstdate='.true.' ; ln_shlat2d='.true.'   ; nn_diaharm=1981
rn_Cd0=2.5e-3       ; ln_loglayer='.true.'  ; ln_tide='.true.'
ln_boost='.false.'  ; ln_gls='.true.'       ; ln_int_wave_drag='.true.'
ln_hpg_djc='.true.' ;                         ln_calc_tdiss='.false.'
########################################################################

#################### runscript options #################################
# Controls for the simulation                                          #
########################################################################
ln_spin='.true.'   # Do we require a spin up phase for the ice
nn_spin=3          # How many years spin up
year_st=1976       # Start year of the simulation (ramains constant)
year_en=1981       # End year of this job submission
offset=0           # If rdt has changed, an offset to nn_it000 is needed
nn_spin_cycle=1    # Change to nn_spin+1 if year is .gt. year_st
year=$year_st      # Change if starting job part way through simulation
########################################################################


########################################################################
#                No need to edit below this line                       #
########################################################################

   #################### nemo runscript options #########################
   # The following options are for the MES domain_cfg.nc               #
   #####################################################################

SRUN_CMD='mpiexec.hydra -print-rank-map -ppn 1 -np 16 ./xios_server.exe : -np 1326 ./nemo'


if [ "$ln_zps" = ".true." ]
then
   ln_trabbl='.true.'
   ln_traldf_hor='.false.'
   ln_traldf_iso='.true.'
   ln_hpg_prj='.false.'
   ln_hpg_sco='.true.'
else
   ln_trabbl='.false.'
   ln_traldf_hor='.true.'
   ln_traldf_iso='.false.'
   if [ "$ln_hpg_djc" = ".true." ]
   then
      ln_hpg_prj='.false.'
   else
      ln_hpg_prj='.true.'
   fi
   ln_hpg_sco='.false.'
fi

if [ "$ln_loglayer" = ".true." ]
then
   ln_non_lin='.false.'
else
   ln_non_lin='.true.'
fi

if [ "$ln_gls" = ".true." ]
then
   ln_tke='.false.'
   ln_zdfevd='.false.'
else
   ln_tke='.true.'
   ln_zdfevd='.true.'
fi

if [ "$ln_rnf_new" = ".true." ]
then
	sn_rnf="'.\/JRA_RIVERS\/ORCA025_rivers_Antar_Green',24,'rorunoff',.false.,.false.,'yearly'"
else
	sn_rnf="'runoff_1m_nomask.nc',-1,'sornficb',.true.,.true.,'yearly'"
fi

echo $SLURM_NNODES

. XXX_ENV_XXX
export OMP_NUM_THREADS=1


jpni=-1
jpnj=-1

nn_date0=$year_st\0101
while [ $year -le $year_en ]
do

   echo $year

   # Array pretending to be a Pythonic dictionary {EXP_NAME:NZ}
   ARRAY=( "SENEMO:75"
         )

   for ens in "${ARRAY[@]}"
   do

      nam="${ens%%:*}"
      jpk="${ens##*:}"

      if [ $year -eq $nn_diaharm ] && [ "$ln_tide" = ".true." ]
      then
	 ln_diaharm='.true.'
      else
	 ln_diaharm='.false.'
      fi

      if [ $year -eq $year_st ]
      then
         nn_it000=1
	 if [ $nn_spin_cycle -eq 1 ]
         then
            ln_rstart='.false.'	  
	    ln_tide_ramp='.true.'
	 else
            ln_rstart='.true.'	  
	    ln_tide_ramp='.false.'
         fi
	 nn_rstctl=0
         ln_reset_ts='.true.'
      else
         nn_it000=`./end_time_step $(( $year - 1 )) $rn_rdt $year_st` ; rs0=`printf "%08d" $nn_it000`; nn_it000=$(( $nn_it000 + 1 + $offset ))
         ln_rstart=".true."
	 nn_rstctl=2
	 ln_tide_ramp='.false.'
         ln_reset_ts='.false.'
      fi
       
      if [ "$ln_rstdate" = ".true." ]
      then
         if [ $year -eq $year_st ]
         then
            rs0=$(($year+1))\0101
         else
            rs0=$year\0101
         fi
      fi 

      nn_itend=`./end_time_step $year $rn_rdt $year_st`; nn_itend=$(( $nn_itend + $offset ))
      if [ $nn_spin_cycle -le $nn_spin ]
      then
         sed    "s/XXX_EXP_XXX/$nam\_S$nn_spin_cycle/g" namelist_cfg_template > namelist_cfg
         sed -i "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle-1))_$rs0\_restart/g" namelist_cfg
         sed "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle-1))_$rs0\_restart_ice/g" namelist_ice_cfg_template > namelist_ice_cfg
      else
         sed    "s/XXX_EXP_XXX/$nam/g" namelist_cfg_template > namelist_cfg
	 if [ $year -eq $year_st ]
         then
            sed -i "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle-1))_$rs0\_restart/g" namelist_cfg
            sed "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle-1))_$rs0\_restart_ice/g" namelist_ice_cfg_template > namelist_ice_cfg
         else
            sed -i "s/XXX_RS0_XXX/$nam\_$rs0\_restart/g" namelist_cfg
            sed "s/XXX_RS0_XXX/$nam\_$rs0\_restart_ice/g" namelist_ice_cfg_template > namelist_ice_cfg
         fi
      fi
      sed -i "s/XXX_TST_XXX/$nn_it000/g" namelist_cfg
      sed -i "s/XXX_YEAR_XXX/$year/g" namelist_cfg
      sed -i "s/XXX_TEN_XXX/$nn_itend/g" namelist_cfg
      sed -i "s/XXX_RDT_XXX/$rn_rdt/g" namelist_cfg
      sed -i "s/XXX_RST_XXX/$ln_rstart/g" namelist_cfg
      sed -i "s/XXX_RCT_XXX/$nn_rstctl/g" namelist_cfg
      sed -i "s/XXX_RSD_XXX/$ln_rstdate/g" namelist_cfg
      sed -i "s/XXX_RTS_XXX/$ln_reset_ts/g" namelist_cfg
      sed -i "s/XXX_RMP_XXX/$ln_tide_ramp/g" namelist_cfg
      sed -i "s/XXX_BTA_XXX/$ln_bt_auto/g" namelist_cfg
      sed -i "s/XXX_BTC_XXX/$rn_bt_cmax/g" namelist_cfg
      sed -i "s/XXX_RNF_XXX/$sn_rnf/g" namelist_cfg
      sed -i "s/XXX_TID_XXX/$ln_tide/g" namelist_cfg
      sed -i "s/XXX_ITF_XXX/$ln_tmx_itf/g" namelist_cfg
      sed -i "s/XXX_DHM_XXX/$ln_diaharm/g" namelist_cfg
      sed -i "s/XXX_BOO_XXX/$ln_boost/g" namelist_cfg
      sed -i "s/XXX_SHA_XXX/$ln_shlat2d/g" namelist_cfg
      sed -i "s/XXX_CD0_XXX/$rn_Cd0/g" namelist_cfg
      sed -i "s/XXX_GLS_XXX/$ln_gls/g" namelist_cfg
      sed -i "s/XXX_TKE_XXX/$ln_tke/g" namelist_cfg
      sed -i "s/XXX_EVD_XXX/$ln_zdfevd/g" namelist_cfg
      sed -i "s/XXX_BBL_XXX/$ln_trabbl/g" namelist_cfg
      sed -i "s/XXX_HOR_XXX/$ln_traldf_hor/g" namelist_cfg
      sed -i "s/XXX_ISO_XXX/$ln_traldf_iso/g" namelist_cfg
      sed -i "s/XXX_PRJ_XXX/$ln_hpg_prj/g" namelist_cfg
      sed -i "s/XXX_DJC_XXX/$ln_hpg_djc/g" namelist_cfg
      sed -i "s/XXX_SCO_XXX/$ln_hpg_sco/g" namelist_cfg
      sed -i "s/XXX_LOG_XXX/$ln_loglayer/g" namelist_cfg
      sed -i "s/XXX_NLN_XXX/$ln_non_lin/g" namelist_cfg
      sed -i "s/XXX_Z0I_XXX/$nn_z0_ice/g" namelist_cfg
      sed -i "s/XXX_IMX_XXX/$nn_mxlice/g" namelist_cfg
      sed -i "s/XXX_PNI_XXX/$nn_jpni/g" namelist_cfg
      sed -i "s/XXX_PNJ_XXX/$nn_jpnj/g" namelist_cfg
      sed -i "s/XXX_TDG_XXX/$ln_int_wave_drag/g" namelist_cfg
      sed -i "s/XXX_CTD_XXX/$ln_calc_tdiss/g" namelist_cfg

   done

   for ens in "${ARRAY[@]}"
   do
      eval $SRUN_CMD &
      echo $SRUN_CMD &
      ./time_step_chk $SLURM_JOB_ID $nn_itend &
   done

   wait
   date

   for ens in "${ARRAY[@]}"
   do
      if [ $year -eq $year_st ] && [ $nn_spin_cycle -le $nn_spin ]
      then
         suf=$year\_S$nn_spin_cycle
      else
         suf=$year
      fi
      #xp="${ens%%:*}"
      #cd $PBS_O_WORKDIR/../ENSEMBLE_MEMBERS/ENS_$xp
      mv *_??_*grid*.nc OUTPUTS &
      mv *_??_*icemod*.nc OUTPUTS &
      cp ocean.output meta_out/ocean.output.$suf
      cp namelist_cfg meta_out/namelist_cfg.$suf
      cp namelist_ice_cfg meta_out/namelist_ice_cfg.$suf
      cp run.stat meta_out/run.stat.$suf
      cp time.step meta_out/time.step.$suf
   done
   wait

   current_stp=`sed -n 1,1p time.step`
   if [ ! $current_stp -eq $nn_itend ]
   then
      exit
   fi

   if [ $nn_spin_cycle -gt $nn_spin ]
   then
      year=$(($year+1))
   else
      nn_spin_cycle=$(($nn_spin_cycle+1))
   fi
      
done

exit
