vcf_dir=/Users/pmonsieurs/programming/leishmania_scDNA_atrandi/results/drugs/

for sample in 1 2 3 5 6; do

    cd /Users/pmonsieurs/programming/software/snpEff
    cd data/

    ## run snpEff
    ln -s ${vcf_dir}/sample_${sample}.filtered.drugs.vcf.gz ${vcf_dir}/sample_${sample}.filtered.drugs.vcf
    java -jar ~/programming/software/snpEff/snpEff.jar Leishmania_donovani_16Nov2015beta  ${vcf_dir}/sample_${sample}.filtered.drugs.vcf > ${vcf_dir}/sample_${sample}.filtered.drugs.snpeff.vcf
    
    ## run snpSift to extract useful information
    cd /Users/pmonsieurs/programming/software/snpEff
    java -jar SnpSift.jar extractFields ${vcf_dir}/sample_${sample}.filtered.drugs.snpeff.vcf CHROM POS "ANN[0].IMPACT" "ANN[0].EFFECT" REF "ANN[0].ALLELE" "ANN[0].HGVS_P" > ${vcf_dir}/sample_${sample}.filtered.drugs.snpeff.csv

done

