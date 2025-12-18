import os
for run_number in range(12):
    runstring=f"""#!/bin/bash
#SBATCH --partition=standard
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --account=class
#SBATCH --qos=standard
cd /home/users/jholt/Git/SE-NEMO/scripts/analysis/Python/

conda activate coast_dev2

python  Physical_indicators_LME.py {run_number}
    """
    with open (f"runscript_inds{run_number}.slurm",'w+') as f:
        f.writelines(runstring)

    os.system(f"sbatch runscript_inds{run_number}.slurm")

