#!/bin/bash
## sbatch -J Thailand lotus_NEMO_profiles_by_LME.sh 34
#SBATCH --partition=short-serial
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=120:00
#SBATCH --ntasks=1
module add jaspy
source activate senemo-profile-lme

echo $1
python  /gws/nopw/j04/class_vol2/senemo/jelt/PROCESSED/EN4.2.1/1978-2019/NEMO_profiles_by_LME.py $1 > LOGS/OUT_$1
