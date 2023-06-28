#!/bin/bash
#SBATCH --job-name=se-025ERSEM
#SBATCH --output=out-se025-%j.out
#SBATCH --time=04:30:00
##SBATCH --time=16:00:00
##SBATCH --time=00:15:00
#SBATCH --account=n01-pml
#SBATCH --partition=standard
#SBATCH --qos=standard
##SBATCH --nodes=19
##SBATCH --nodes=52
#SBATCH --nodes=68
##SBATCH --nodes=99
#SBATCH --ntasks-per-core=1

# ARCHER2 support suggestion to reduce UCX error messages:
export UCX_IB_REG_METHODS=direct

RUNDIR=$PWD

# SUBMIT WITH:
# "sbatch --export=y0=1976,m0=1,nn_spin_cycle=0 runscript.multimonth.slurm"
# NB: m0 is month to run first, 
#     nn_spin_cycle is spin cycle of the previous month (i.e. that tagged in the restart)
#     After the spin cycle it should be sufficient to have any nn_spin_cycle>=nn_spin+1 (>=4)


#################### nemo runscript options ############################
# For info on the parameters see namelist_ref                          #
########################################################################
rn_rdt=1800          ; ln_zps='.true.'       ; ln_tmx_itf='.true.'
ln_bt_auto='.true.' ; rn_bt_cmax=0.8        ; nn_baro=30
nn_mxlice=3         ; nn_z0_ice=1           ; ln_rnf_new='.true.'
ln_rstdate='.false.' ; ln_shlat2d='.false.'  ; nn_diaharm=3333
rn_Cd0=1.0e-3       ; ln_loglayer='.false.' ; ln_tide='.false.'
ln_boost='.true.'   ; ln_gls='.false.'       ; ln_int_wave_drag='.false.' 
ln_hpg_djc='.false.' ;
########################################################################


#################### runscript options #################################
# Controls for the simulation                                          #
########################################################################
ln_spin='.true.'                # Do we require a spin up phase for the ice
nn_spin=3                       # How many years spin up
year_st=1976                    # Start year of the simulation (ramains constant)
year_en=2019                    # End year of this job submission
offset=0                        # If rdt has changed, an offset to nn_it000 is needed
nn_spin_cycle=$nn_spin_cycle    # Change to nn_spin+1 if year is .gt. year_st
year=$y0                        # Change if starting job part way through simulation
nmmax=12                        # number of months to go through before submitting a new job

########################################################################


########################################################################
#                No need to edit below this line                       #
########################################################################
# this has been added because on the new opearting system quite often if the file exists it is not overwritten
if [ -f "myscript_wrapper.sh" ]; then
        rm myscript_wrapper.sh
fi

# this has been added to ensure that right domain file is linked
if [ -f "domain_cfg.nc" ]; then
        rm domain_cfg.nc
fi
if [ "$ln_zps" = ".true." ]
then
        ln -s /work/n01/shared/yuti/se-ORCA025_domain_zps/domain_cfg_zps.closea.nc domain_cfg.nc
else
        ln -s INPUTS/domain_cfg_mes_v2.nc ./domain_cfg.nc
fi

if [ "$ln_zps" = ".false." ]
then

   #################### nemo runscript options #########################
   # The following options are for the MES domain_cfg.nc               #
   #####################################################################

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
   elif [ $SLURM_NNODES -eq 99 ]
   then
      SRUN_CMD=`./slurm_setup -S 32 -s 16 -m 2 -C 9104 -g 4 -N 128 -t 01:00:00 -a n01-CLASS -j se-nemo`
   else
      exit
   fi

else

   #################### nemo runscript options #########################
   # The following options are for the ZPS domain_cfg.nc               #
   #####################################################################

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
   elif [ $SLURM_NNODES -eq 99 ]
   then
      SRUN_CMD=`./slurm_setup -S 32 -s 16 -m 2 -C 9104 -g 4 -N 128 -t 01:00:00 -a n01-CLASS -j se-nemo`
   else
      exit
   fi

fi


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
	cn_dir=""
	sn_rnf="'JRA_BGC',24,'rorunoff',.true.,.false.,'yearly'"
	sn_cnf="'JRA_BGC',0,'socoefr',.false.,.false.,'yearly'"
	sn_dep="'JRA_BGC',24,'rodepth',.false.,.false.,'yearly'"
else
	cn_dir="INPUTS"
	sn_rnf="'runoff_1m_nomask.nc',-1,'sornficb',.true.,.true.,'yearly'"
	sn_cnf="'runoff_1m_nomask.nc',0,'socoefr',.false.,.true.,'yearly'"
fi

echo $SLURM_NNODES

. /work/n01/n01/yuti/SE-NEMO/scripts/env/gnu-mpich
export OMP_NUM_THREADS=1

chmod u+x ./myscript_wrapper.sh

jpni=-1
jpnj=-1

nsstart=$(date -d "$year_st-01-01 00:00:00 GMT" +%s) #seconds since EPOCH for total simulation start


#set annual links for ERSEM
./set_annual_links.sh $y0

stopflag=0

nm=1
while [ $y0 -le $year_en ] && [ $nm -le $nmmax ]
do

   rm timing_$nm.txt
   touch timing_$nm.txt
   echo start: >> timing_$nm.txt
   ls -lhrt timing_$nm.txt >> timing_$nm.txt

   #Compute previous and next month:
   mp=$(( $m0 + 1 ))
   if [ $mp -eq 13 ]
   then
      mp=1
      if [ $nn_spin_cycle -gt $nn_spin ]
      then
        yp=$(( y0 + 1 ))
      else
	yp=$y0
      fi
   else
      yp=$y0
   fi
   mm=$(( $m0 - 1 ))
   if [ $mm -eq 0 ]
   then
      mm=12
      if [ $nn_spin_cycle -gt $nn_spin ]
      then
          ym=$(( y0 - 1 ))
      else
	  ym=$y0
	  nn_spin_cycle=$(($nn_spin_cycle+1))
      fi
   else
      ym=$y0
   fi
   echo "Preparing $y0 $m0 ..."
   #./monthlyclean.sh $RUNDIR   #TO DO!!!!!!!!!!!!!!

   m0str=$(printf %02d $m0)
   mpstr=$(printf %02d $mp)
   mmstr=$(printf %02d $mm)

   # Array pretending to be a Pythonic dictionary {EXP_NAME:NZ}
   ARRAY=( "SENEMO:75"
         )

   for ens in "${ARRAY[@]}"
   do

      nam="${ens%%:*}"
      jpk="${ens##*:}"

      if [ $y0 -eq $nn_diaharm ]
      then
	 ln_diaharm='.true.'
      else
	 ln_diaharm='.false.'
      fi

      ln_rstdate='.false.'
      if [ $y0 -eq 1976 ]
      then
	if [ $m0 -eq 1 ]
	then
	  ln_rstdate='.false.'
	fi
      fi

      if [ $y0 -eq 1976 ]
      then
         nn_it000=1
	 if [ $nn_spin_cycle -eq 1 ]
         then
            ln_rstart='.false.'	  
            ln_rstart_top='.true.'	  
	    ln_tide_ramp='.true.'
	 else
            ln_rstart='.true.'	  
            ln_rstart_top='.true.'
	    ln_tide_ramp='.false.'
         fi
	 nn_rstctl=0
         ln_reset_ts='.true.'
      else
         nn_it000=`./end_time_step $(( $y0 - 1 )) $rn_rdt $year_st` ; rs0=`printf "%08d" $nn_it000`; nn_it000=$(( $nn_it000 + 1 + $offset ))
         ln_rstart=".true."
         ln_rstart_top=".true."
	 nn_rstctl=2
	 ln_tide_ramp='.false.'
         ln_reset_ts='.false.'
      fi

       # number of days in the current month
       case $m0 in
      	       4|6|9|11) ndays=30 ;;
      	       2) if [ $(( y0 % 4 )) -ne 0 -o $(( y0 % 100)) -eq 0 -a $(( $y0 % 400 )) -ne 0 ]; then ndays=28; else ndays=29; fi ;;
      	       *) ndays=31 ;;
       esac
       # number of days in the previous month (needed for restart_date)
       case $mm in
      	       4|6|9|11) mdays=30 ;;
      	       2) if [ $(( ym % 4 )) -ne 0 -o $(( ym % 100)) -eq 0 -a $(( $ym % 400 )) -ne 0 ]; then mdays=28; else mdays=29; fi ;;
      	       *) mdays=31 ;;
       esac

       echo 'CALCULATE NIT'

       nit=$(( 86400 * $ndays / $rn_rdt ))

       echo 'nit= '$nit

       ns0=$(date -d "$y0-${m0str}-01 00:00:00 GMT" +%s) #seconds since EPOCH for this chunk
       nn_it000=$(( $ns0 - $nsstart ))
       nn_it000=$(( $nn_it000 / $rn_rdt + 1 + $offset ))
       nn_itend=$(( $nn_it000 + $nit -1 ))
       nn_date0=$y0${m0str}01
       date_restart=$ym${mmstr}$mdays
       echo nn_it000: $nn_it000 !!!!!

      nn_itrst=$(( nn_it000 - 1 ))
      rs0=$(printf "%08d" $nn_itrst)
      echo rs0: $rs0 !!!!

      nn_spin_cycle_rst=$nn_spin_cycle
      if [ $m0 -eq 1 ]
      then
        nn_spin_cycle_rst=$(( nn_spin_cycle - 1 ))
      fi

      #set physical namelist
      if [ $nn_spin_cycle -le $nn_spin ]
      then
         sed    "s/XXX_EXP_XXX/$nam\_S$nn_spin_cycle/g" namelist_cfg_template > namelist_cfg   
         sed -i "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle_rst))_$rs0\_restart/g" namelist_cfg    # restart file NAME of the member_SX_YYYY_restart
         sed "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle_rst))_$rs0\_restart_ice/g" namelist_ice_cfg_template > namelist_ice_cfg
         sed "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle_rst))_$rs0\_restart_trc/g" namelist_top_template > namelist_top_cfg
      else
         sed    "s/XXX_EXP_XXX/$nam/g" namelist_cfg_template > namelist_cfg
         if [ $y0 -eq $year_st ] && [ $m0 -eq 1 ]
         then
            sed -i "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle_rst))_$rs0\_restart/g" namelist_cfg  
            sed "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle_rst))_$rs0\_restart_ice/g" namelist_ice_cfg_template > namelist_ice_cfg
            sed "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle_rst))_$rs0\_restart_trc/g" namelist_top_template > namelist_top_cfg
         else
            sed -i "s/XXX_RS0_XXX/$nam\_$rs0\_restart/g" namelist_cfg
            sed "s/XXX_RS0_XXX/$nam\_$rs0\_restart_ice/g" namelist_ice_cfg_template > namelist_ice_cfg
            sed "s/XXX_RS0_XXX/$nam\_$rs0\_restart_trc/g" namelist_top_template > namelist_top_cfg
         fi
      fi

      sed -i "s/XXX_TST_XXX/$nn_it000/g" namelist_cfg
      #sed -i "s/XXX_YEAR_XXX/$year/g" namelist_cfg
      sed -i "s/XXX_TEN_XXX/$nn_itend/g" namelist_cfg
      sed -i "s/XXX_RDT_XXX/$rn_rdt/g" namelist_cfg
      sed -i "s/XXX_RST_XXX/$ln_rstart/g" namelist_cfg
      sed -i "s/XXX_RCT_XXX/$nn_rstctl/g" namelist_cfg
      sed -i "s/XXX_RSD_XXX/$ln_rstdate/g" namelist_cfg
      sed -i "s/XXX_RTS_XXX/$ln_reset_ts/g" namelist_cfg
      sed -i "s/XXX_DATE0_XXX/$nn_date0/g" namelist_cfg
      sed -i "s/XXX_RMP_XXX/$ln_tide_ramp/g" namelist_cfg
      sed -i "s/XXX_BTA_XXX/$ln_bt_auto/g" namelist_cfg
      sed -i "s/XXX_BTC_XXX/$rn_bt_cmax/g" namelist_cfg
      sed -i "s/XXX_CNDIR_XXX/$cn_dir/g" namelist_cfg
      sed -i "s/XXX_RNF_XXX/$sn_rnf/g" namelist_cfg
      sed -i "s/XXX_CNF_XXX/$sn_cnf/g" namelist_cfg
      sed -i "s/XXX_DEP_XXX/$sn_dep/g" namelist_cfg
      sed -i "s/XXX_TID_XXX/$ln_tide/g" namelist_cfg
      sed -i "s/XXX_ITF_XXX/$ln_tmx_itf/g" namelist_cfg
      sed -i "s/XXX_DHM_XXX/$ln_diaharm/g" namelist_cfg
      sed -i "s/XXX_BOO_XXX/$ln_boost/g" namelist_cfg
      sed -i "s/XXX_DJC_XXX/$ln_hpg_djc/g" namelist_cfg
      sed -i "s/XXX_SHA_XXX/$ln_shlat2d/g" namelist_cfg
      sed -i "s/XXX_CD0_XXX/$rn_Cd0/g" namelist_cfg
      sed -i "s/XXX_GLS_XXX/$ln_gls/g" namelist_cfg
      sed -i "s/XXX_TKE_XXX/$ln_tke/g" namelist_cfg
      sed -i "s/XXX_EVD_XXX/$ln_zdfevd/g" namelist_cfg
      sed -i "s/XXX_BBL_XXX/$ln_trabbl/g" namelist_cfg
      sed -i "s/XXX_HOR_XXX/$ln_traldf_hor/g" namelist_cfg
      sed -i "s/XXX_ISO_XXX/$ln_traldf_iso/g" namelist_cfg
      sed -i "s/XXX_PRJ_XXX/$ln_hpg_prj/g" namelist_cfg
      sed -i "s/XXX_SCO_XXX/$ln_hpg_sco/g" namelist_cfg
      sed -i "s/XXX_LOG_XXX/$ln_loglayer/g" namelist_cfg
      sed -i "s/XXX_NLN_XXX/$ln_non_lin/g" namelist_cfg
      sed -i "s/XXX_Z0I_XXX/$nn_z0_ice/g" namelist_cfg
      sed -i "s/XXX_IMX_XXX/$nn_mxlice/g" namelist_cfg
      sed -i "s/XXX_PNI_XXX/$nn_jpni/g" namelist_cfg
      sed -i "s/XXX_PNJ_XXX/$nn_jpnj/g" namelist_cfg
      sed -i "s/XXX_TDG_XXX/$ln_int_wave_drag/g" namelist_cfg
      
      sed -i "s/XXX_RST_XXX/$ln_rstart_top/g" namelist_top_cfg
      sed -i "s/XXX_RCT_XXX/$nn_rstctl/g" namelist_top_cfg

   done 

   # execution
   for ens in "${ARRAY[@]}"
   do
      ##time_step_chk sometimes interferes with execution, makes the run hang.
      #eval $SRUN_CMD &
      #echo $SRUN_CMD &
      #./time_step_chk $SLURM_JOB_ID $nn_itend &
      eval $SRUN_CMD &
      echo $SRUN_CMD 
   done

   wait
   echo "Finished $y0 $m0 at $(date +%s) seconds since 1970-01-01 00:00:00"

   for ens in "${ARRAY[@]}"
   do
      if [ $y0 -eq $year_st ] && [ $nn_spin_cycle -le $nn_spin ]
      then
         suf=$y0$m0\_S$nn_spin_cycle
	 res=$y0$m0$ndays\_S$nn_spin_cycle
      else
         suf=$y0$m0
	 res=$y0$m0$ndays
      fi

      #prepare archive directory
      ARCHIVEDIR=$RUNDIR/$nam/$y0/$m0str/
      mkdir -p $ARCHIVEDIR

      # compress and archive output
      sbatch --export=ARCHIVEDIR=$ARCHIVEDIR,y0=$y0,m0str=$m0str,SPIN=$nn_spin_cycle compress_output.slurm

      mkdir -p meta_out
      mv -f ocean.output meta_out/ocean.output.bz2.$suf
      cp namelist_cfg meta_out/namelist_cfg.$suf
      cp namelist_ice_cfg meta_out/namelist_ice_cfg.$suf
      cp namelist_top_cfg meta_out/namelist_top_cfg.$suf
      cp run.stat meta_out/run.stat.$suf
      cp time.step meta_out/time.step.$suf

      # REBUILD and COMPRESS RESTART from previous 2 months
      nit_mm1=$(( 86400 * $mdays / $rn_rdt ))
      nn_itend_mm2=$(( $nn_it000 - $nit_mm1 -1 ))
      echo nn_itend_mm2: $nn_itend_mm2

      # compress and rebuild 2 restarts per year (December, June), else delete
      nn_spin_cycle_rbld_rst=$nn_spin_cycle
      if [ $m0 -le 2 ]
      then
        nn_spin_cycle_rbld_rst=$(( nn_spin_cycle_rbld_rst - 1 ))
      fi
      if [ $m0 -eq 2 ] || [ $m0 -eq 8 ]
      then
        rbld_id=$(sbatch --parsable --export=ITER=$nn_itend_mm2,SPIN=$nn_spin_cycle_rbld_rst rebuild_restart.slurm)
        sbatch --dependency=afterok:$rbld_id --export=ITER=$nn_itend_mm2,SPIN=$nn_spin_cycle_rbld_rst compress_restarts.slurm
      else
        bash clean_restarts.sh $nn_itend_mm2 $nn_spin_cycle_rbld_rst
      fi


      # during the spin cycle rename last restart of each year so it sees it
      echo current parameters:
      echo m0: $m0, rs0: $rs0, nn_spin_cycle: $nn_spin_cycle, nn_spin: $nn_spin
      if [ $m0 -eq 12 ] && [ $nn_spin_cycle -le $nn_spin ]
      then
	  echo renaming December restarts
          rename $rs0 00000000 RESTARTS/${nam}_S${nn_spin_cycle}_${rs0}_restart*_????.nc
      fi
   done
   wait

   current_stp=`sed -n 1,1p time.step`
   if [ ! $current_stp -eq $nn_itend ]
   then
      exit
   fi

   echo done 1 month: >> timing_$nm.txt
   touch timing_$nm.txt
   ls -lhrt timing_$nm.txt >> timing_$nm.txt

   y0=$yp
   m0=$mp
   nm=$(( nm + 1 ))

done



##resubmit:
if [ $y0 -le $year_en ]
then
   echo 'Submitting' $y0 $m0 '...'
   sbatch --export=y0=$y0,m0=$m0,nn_spin_cycle=$nn_spin_cycle runscript.multimonth.slurm
   echo "Done."
else
   echo "All done."
fi

exit
