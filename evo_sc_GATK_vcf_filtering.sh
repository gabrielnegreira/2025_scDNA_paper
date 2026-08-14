module load BioTools

cd /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/evo/gvcf/

## select only SNPs with a minor allele frequency (MAF) of at least 0.05 and a missing genotype rate of at most 0.05
for strain in "HU3" "BPK081"; do
    for sample in 1 2 3 5 6; do
        bcftools view -i 'MAF[0] >= 0.05 && F_MISSING <= 0.25' ${strain}/sample_${sample}.filtered.snps.${strain}.vcf.gz -Oz -o ${strain}/sample_${sample}.filtered.snps.${strain}.clean.vcf.gz
    done
done
