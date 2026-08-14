#!/usr/bin/env python3

# module load Python/3.8.3-intel-2020a
 
import pandas as pd
import os

src_dir = '/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/bwa/'
out_dir = '/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/evo/gvcf/'


meta_data_file = '/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/data/cells_meta.xlsx'
meta_data = pd.read_excel(meta_data_file)
print(meta_data)

for index, row in meta_data.iterrows():

    ## extract all the required information
    batch = row['library']
    barcode = f"{row['barcode_D']}_{row['barcode_C']}_{row['barcode_B']}_{row['barcode_A']}"
    sample = row['sample']
    

    ## link the g.vcf file to the corresponding directory
    gvcf_file_in = f"{src_dir}/{batch}/{barcode}.gatk.g.vcf"
    gvcf_file_out = f"{out_dir}/sample_{sample}/{barcode}.gatk.g.vcf"
    ln_command = f"ln -s {gvcf_file_in} {gvcf_file_out}"
    print(ln_command)
    os.system(ln_command)