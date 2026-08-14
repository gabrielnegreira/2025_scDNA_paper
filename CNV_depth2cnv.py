#!/usr/bin/env python3

import pandas as pd
import argparse


## the lines in the initial gff file are not the intended number of lines, so this
## should be removed. Those are the Ld37 genes that indeed contain too much information
## on those lines
# awk -F'\t' 'NF == 9' TriTrypDB-64_LtropicaL590.gff
# grep -v "#"  TriTrypDB-64_LtropicaL590.gff > TriTrypDB-64_LtropicaL590.clean.gff
# grep "exon"  TriTrypDB-64_LtropicaL590.clean.gff > TriTrypDB-64_LtropicaL590.clean.exon.gff

data_dir = '/user/antwerpen/205/vsc20587/scratch/leishmania_susl/data/refgenomes/'
gff_file = f"{data_dir}/Leishmania_donovani_16Nov2015beta.clean.gff3"

## function to read in the gff data
def read_gff (gff_file):
    print(gff_file)
    try:
        gff_data = pd.read_csv(gff_file, sep='\t', header=None)
        gff_data = gff_data[[0, 3, 4, 8]]
        gff_data[8] = gff_data[8].str.extract(r'(LdBPK_\d+)')
        gff_data.columns = ['chrom', 'start', 'end', 'gene_id']
        print("GFF Data Loaded Successfully")
        return gff_data
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        return None

def read_coverage (cov_file):
    print(cov_file)
    try:
        cov_data = pd.read_csv(cov_file, sep='\t', header=None)
        cov_data.columns = ['chrom', 'pos', 'cov']
        print("Coverage Data Loaded Successfully")
        return cov_data
    except Exception as e:
        print(f"Error reading Coverage file: {e}")
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculate the relative coverage of a gene based on the output of samtools depth command (.cov file) combined with the .gff file defining the genes")
    parser.add_argument("-c", "--cov_file", help="path to the file containing the output of the samtools depth command")

    ## parse arguments
    args = parser.parse_args()
    cov_file = args.cov_file
    print(cov_file)

    ## parse the gff data
    gff_data = read_gff(gff_file)

    ## parse the coverage data 
    print(cov_file)
    cov_data = read_coverage(cov_file)
    median_coverage_per_chromosome = cov_data.groupby("chrom")['cov'].median()
    print(median_coverage_per_chromosome)


    gene_cnv = []

    count = 0 
    for _, gene in gff_data.iterrows():

        ## extract relevant information from gff_data
        chrom = gene['chrom']
        start = gene['start']
        end = gene['end']

        print(f"{chrom}:{start}-{end}")
        ## filter coverage data to get the relevant rows for the gene (chromosome and positions within start and end) and get the median value
        gene_cov_data = cov_data[(cov_data['chrom'] == chrom) & (cov_data['pos'] >= start) & (cov_data['pos'] <= end)]
        gene_median_cov = gene_cov_data['cov'].median()
        
        ## retrieve the corresponding median for the chromosome, and do the ratio
        median_cov_chrom = median_coverage_per_chromosome[chrom]
        ratio = gene_median_cov / median_cov_chrom if median_cov_chrom != 0 else "NA"  # Prevent division by zero

        ## store the results
        gene_cnv.append({'gene_id': gene['gene_id'], 'median_coverage': gene_median_cov, 'chrom': chrom, 'chrom_median': median_cov_chrom, 'ratio': ratio})

        ## do some printing to check
        print(f"gene_id {gene['gene_id']} - median coverage {gene_median_cov} - chrom_median: {median_cov_chrom} - ratio {ratio}")

    ## convert the results to a DataFrame and save to .csv
    gene_cnv_df = pd.DataFrame(gene_cnv)
    print(gene_cnv_df)
    
    output_file = cov_file.replace(".depth", ".cnv.csv")
    gene_cnv_df.to_csv(output_file, index=False)





