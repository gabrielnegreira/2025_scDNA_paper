#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=50:00:00
#SBATCH --job-name=scistree2_run
#SBATCH --account=ap_itg_mpu
#SBATCH --mail-user=gnegreira@itg.be
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --array=1-NUM_JOBS  # will be replaced dynamically

# Load modules
module --force purge
module load calcua/all
module load hpc-container-wrapper

# Activate container
export PATH="$VSC_SCRATCH/containers/scistree2/bin:$PATH"
unset PYTHONPATH

# Create a sorted list of files
FILE_LIST=($(ls "$VSC_SCRATCH/atrandi_paper/inputs/nucleotide_variants"/*for_scistree2.tsv))

# Select the file corresponding to this array task
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing file: $FILE"
python ScisTree2.py "$FILE"