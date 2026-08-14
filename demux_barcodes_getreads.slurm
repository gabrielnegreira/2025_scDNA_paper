#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=50:50:00
#SBATCH --job-name=demux_fastq
#SBATCH -A ap_itg_mpu

module load atools
module load Python/3.12.3-GCCcore-13.3.0

## export PATH and specify parameters 
export PYTHONPATH="/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.12/site-packages/:${PYTHONPATH}"
export PATH=/user/antwerpen/205/vsc20587/data/software/python_lib/bin/:/user/antwerpen/205/vsc20587/data/software/seqtk/:/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/:${PATH}

#### scDNA_AT_01 ####
# demux_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/scDNA_AT_01/
# fastq_file_R1=/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/data/fastq/scDNA_AT_01_1.fq.gz
# fastq_file_R2=/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/data/fastq/scDNA_AT_01_2.fq.gz
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/scDNA_AT_01_barcodes.txt)

#### scDNA_AT_02 ####
demux_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/scDNA_AT_02/
fastq_file_R1=/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/data/fastq/scDNA_AT_02_1.fq.gz
fastq_file_R2=/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/data/fastq/scDNA_AT_02_2.fq.gz
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/scDNA_AT_02_barcodes.txt)
echo "barcdode ${barcode}"

## extract for the barcode all the read IDs that should be extracted
## from the fastq-files. The output of this 
demux_barcodes_getreads.py --outputdir ${demux_dir} --fastq ${fastq_file_R2} --barcode ${barcode}

## construct the output file name of the demux_barcodes_getreads.py command
## and create the output files of the demultiplexing
read_id_file=${demux_dir}/${barcode}.readids.txt
fastq_file_barcode_R1=${demux_dir}/${barcode}_R1.fastq.gz
fastq_file_barcode_R2=${demux_dir}/${barcode}_R2.fastq
fastq_file_barcode_R2_with_bc=${demux_dir}/${barcode}_R2_with_bc.fastq.gz


seqtk subseq ${fastq_file_R1} ${read_id_file} | gzip > ${fastq_file_barcode_R1}
seqtk subseq ${fastq_file_R2} ${read_id_file} | gzip > ${fastq_file_barcode_R2_with_bc}
cutadapt --cut 48 -o ${fastq_file_barcode_R2} ${fastq_file_barcode_R2_with_bc}
gzip ${fastq_file_barcode_R2}


## run for all barcodes
## sample scDNA_01
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/
# cut -f 1 scDNA_AT_01_barcode_counts.txt > scDNA_AT_01_barcodes.txt
# vi /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/scDNA_AT_01_barcodes.txt ## add barcode to the top of the file
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/scDNA_AT_01_barcodes.txt) /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/demux_barcodes_getreads.slurm


## sample scDNA_02
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/
# cut -f 1 scDNA_AT_02_barcode_counts.txt > scDNA_AT_02_barcodes.txt
# vi /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/scDNA_AT_02_barcodes.txt ## add barcode to the top of the file
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux/scDNA_AT_02_barcodes.txt) /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/demux_barcodes_getreads.slurm



## installation of cutadapt
# module load Python/3.12.3-GCCcore-13.3.0
# cd /user/antwerpen/205/vsc20587/data/software/python_lib/lib
# mkdir -p python3.12/site-packages/
# export PYTHONPATH="/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.12/site-packages/:${PYTHONPATH}"
# pip3 install --prefix=/user/antwerpen/205/vsc20587/data/software/python_lib/ cutadapt


