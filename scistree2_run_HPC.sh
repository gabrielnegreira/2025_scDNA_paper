#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=64
#SBATCH --time=50:00:00
#SBATCH --job-name=scistree2_run
#SBATCH --account=ap_itg_mpu
#SBATCH --mail-user=gnegreira@itg.be
#SBATCH --mail-type=BEGIN,END,FAIL



#load modules
module load calcua/all
module load hpc-container-wrapper

#activate container
export PATH="$VSC_SCRATCH/containers/scistree2/bin:$PATH"
unset PYTHONPATH

#run the script
