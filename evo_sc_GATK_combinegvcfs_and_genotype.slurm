#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=22:30:00
#SBATCH --job-name=bwa_gatk
#SBATCH -A ap_itg_mpu


module load GATK
module load BCFtools
module load Java

## sample should be read from the command line! using --export
ref_genome=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/refgenomes/Leishmania_donovani_16Nov2015beta.fa
vcf_dir=/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/evo/gvcf/sample_${sample}/
gvcf_output=${vcf_dir}/sample_${sample}.g.vcf.gz


## combine gvcf
cd $vcf_dir
gatk CombineGVCFs  -R ${ref_genome}  $(for vcf in *.vcf; do  echo "-V ${vcf} "; done) -O ${gvcf_output}


## create index file and do the genotyping
gatk IndexFeatureFile --input ${vcf_dir}/sample_${sample}.g.vcf.gz
gatk GenotypeGVCFs -R ${ref_genome} -V ${vcf_dir}/sample_${sample}.g.vcf.gz -O ${vcf_dir}/sample_${sample}.vcf.gz


## split up in SNP and indels, and do separate
## filtering at both levels


cd ${vcf_dir}

##### SNPs #####
## select SNPs and do basic filtering on the SNPs based
## on the GATK recommended settings
gatk SelectVariants -V ${vcf_dir}/sample_${sample}.vcf.gz  --select-type-to-include SNP  -O ${vcf_dir}/sample_${sample}_snps.vcf.gz

## add filter information
gatk VariantFiltration -V ${vcf_dir}/sample_${sample}_snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2"  \
    -filter "QUAL < 200.0" --filter-name "QUAL200" \
    -filter "SOR > 3.0" --filter-name "SOR3"  \
    -filter "FS > 60.0" --filter-name "FS60"  \
    -filter "MQ < 40.0" --filter-name "MQ40"  \
    -O ${vcf_dir}/sample_${sample}_filter_added_snps.vcf.gz


## select only SNPs that passed the filter
gatk SelectVariants -R ${ref_genome} -V ${vcf_dir}/sample_${sample}_filter_added_snps.vcf.gz --exclude-filtered true -O ${vcf_dir}/sample_${sample}.filtered.snps.vcf.gz
zcat ${vcf_dir}/sample_${sample}.filtered.snps.vcf.gz | grep -v "#" | wc



##### indels #####
gatk SelectVariants -V ${vcf_dir}/sample_${sample}.vcf.gz -select-type INDEL -O ${vcf_dir}/sample_${sample}.indels.vcf.gz
zcat ${vcf_dir}/${prefix}.indels.vcf.gz | grep -v #" | wc

## add filter information
gatk VariantFiltration  -V ${vcf_dir}/sample_${sample}.indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 200.0" --filter-name "QUAL200" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -O ${vcf_dir}/sample_${sample}.filter_added.indels.vcf.gz

## select only indels that passed the filter - 5347 passed filter
gatk SelectVariants -R $ref_genome -V ${vcf_dir}/sample_${sample}.filter_added.indels.vcf.gz --exclude-filtered true -O ${vcf_dir}/sample_${sample}.filtered.indels.vcf.gz

zcat ${vcf_dir}/sample_${sample}.filtered.indels.vcf.gz | grep -v #" | wc


## merge the SNPs and Indels again to produce one output file
java -jar ~/data/software/picard/picard.jar  SortVcf I=${vcf_dir}/sample_${sample}.filtered.snps.vcf.gz I=${vcf_dir}/sample_${sample}.filtered.indels.vcf.gz O=${vcf_dir}/sample_${sample}.filtered.vcf.gz


## run for different samples
# for sample in 1 2 3 5 6; do
#     sbatch --export=sample=${sample} /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/bin/evo_sc_GATK_combinegvcfs_and_genotype.slurm
# done
