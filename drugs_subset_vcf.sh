module load BioTools

drug_resistance_file=/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/drug_resistance_ldon_genes.bed

cd /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/evo/gvcf/

for sample in 1 2 3 5 6; do
    bcftools view -R ${drug_resistance_file} sample_${sample}/sample_${sample}.filtered.vcf.gz -o sample_${sample}.filtered.drugs.vcf.gz
done







