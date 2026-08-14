#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=30:40:00
#SBATCH -A ap_itg_mpu

module load Python/3.8.3-intel-2020a
module load atools

source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa_depth/depth_files.txt)

/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/CNV_depth2cnv.py -c ${cov_file}

## loop over cov file
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa_depth/
# option1: only taking the filtered cells from Gabriel from cells_meta_v2.xlsx
# vi depth_files.txt
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa_depth/depth_files.txt) /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/CNV_depth2cnv_wrapper.slurm

# option2: looping over all depth files
# for cov_file in *.depth; do sbatch --export=cov_file=${PWD}/${cov_file} /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/CNV_depth2cnv_wrapper.slurm; done

