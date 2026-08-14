#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=7
#SBATCH --time=04:00:00
#SBATCH -A ap_itg_mpu


## installation of deeptools: apparently the -I option was needed to force re-installation of some 
## libraries. Not clear whether this was a deeptools issue, or an issue with a previous 
## installation process
# module load Python/3
# export PYTHONPATH="/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.8/site-packages/:${PYTHONPATH}"
# pip3 install -I  --prefix=/user/antwerpen/205/vsc20587/data/software/python_lib/ deeptools


## specify the deeptools binary and export the PythonPath so it can be run
## with a local binfile
module load Python/3.10.4-GCCcore-11.3.0
module load atools
module load BioTools
export PYTHONPATH="/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.10/site-packages/:${PYTHONPATH}"
bamcoverage_bin=/user/antwerpen/205/vsc20587/data/software/python_lib/bin/bamCoverage

## load the bam-files using the atools option
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/bam_files_sample${sample}_top${top}.txt)
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/bam_files_selection.csv)



## normalization procedure: CPM should work, as this normalizes for the difference
## in sequencing depth. It does not take into account the feature length (in case of RNA-seq)
## but as the main aim is to check the sequencing depth per bin, and the binsize is a constant
## factor, CPM should be sufficient, and is easy interpretable.
## ==> --normalizeUsing CPM

echo $bam_file
samtools index $bam_file

## parameters for running & creating of name for output
threads=7
bin_size=100
output_file=${bam_file/.mapq30.removedups.proper_paired.bam/.deeptools.bin${bin_size}.full_genome.bedgraph}

output_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/cnv/
output_file="${output_file##*/}"
output_file=${output_dir}/all/${output_file}
echo $output_file

$bamcoverage_bin \
    --bam $bam_file \
    --outFileName $output_file \
    --outFileFormat bedgraph \
    --binSize $bin_size \
    --numberOfProcessors $threads \
    --normalizeUsing None
    # --region Ld36:2558141:2576070 \
    # --region Ld36:2300000:2800000 \
    # --region Ld29:961792:963222 \
    # --region Ld29:940000:980000 \


## run for bam files with the highest number of reads per cell - this list is obtained by
## running evo_sc_Mlocus_select_cells.py
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/
# sample=5
# top=154
# sample=6
# top=242
# module load atools
# sbatch --export=sample=${sample},top=${top} --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/bam_files_sample${sample}_top${top}.txt) /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/evo_sc_deeptools.slurm

## option2: run for all the bam_files that have been selected by Gabriel. Create the full path 
## to the .bam file in the excel file "cells_meta_v2.xlsx" and copy paste to new bam_files.csv
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa
# vi bam_files_selection.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/bam_files_selection.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/evo_sc_deeptools.slurm

