#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=70:30:00
#SBATCH -A ap_itg_mpu

# load necessary software settings for GATK
module load Java
module load GATK
module load BioTools
module load atools

## load bam files via sbatch command
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/scDNA_AT_01/bam_files.csv)
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/scDNA_AT_02/bam_files.csv)



## hard coded example
# bam_file=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/results/bwa_plosntds/VL9.mapq30.removedups.proper_paired.bam
ref_genome=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/refgenomes/Leishmania_donovani_16Nov2015beta.fa

## parameter settings
threads=1

## index the bam file - if needed
samtools index ${bam_file}

## GATK using the VCF approach (no population genomics) --> only for diagnostic SNPs
# vcf_file=$(echo $bam_file | sed -e 's/.mapq30.removedups.proper_paired.bam$/.gatk.vcf/')
# gatk HaplotypeCaller -R $ref_genome -I $bam_file -O $vcf_file --native-pair-hmm-threads $threads

## GATK using the GVCF approach (for population genomics approach) 
vcf_file=$(echo $bam_file | sed -e 's/.mapq30.removedups.proper_paired.bam$/.gatk.g.vcf/')
gatk HaplotypeCaller -R $ref_genome -ERC GVCF -I $bam_file -O $vcf_file --native-pair-hmm-threads $threads


## eventually: subset vcf to a list of diagnostic SNPs using bcftools? to be explored


## run for all bam files of batch 1
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/scDNA_AT_01/
# ls $PWD/*paired.bam > bam_files.csv
# vi /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/scDNA_AT_01/bam_files.csv ## add header
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/scDNA_AT_01/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/evo_sc_GATK_haplotypecaller.slurm

## run for all bam files of batch 2
# cd /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/scDNA_AT_02/
# ls $PWD/*paired.bam > bam_files.csv
# vi /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/scDNA_AT_02/bam_files.csv ## add header
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/scDNA_AT_02/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/evo_sc_GATK_haplotypecaller.slurm