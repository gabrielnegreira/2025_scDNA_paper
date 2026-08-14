#!/usr/bin/env python3

import gzip
from collections import defaultdict
import argparse


def hamming_distance(str1, str2):

    if len(str1) != len(str2):
        return False

    # Count the number of differing characters
    differences = sum(1 for a, b in zip(str1, str2) if a != b)

    return differences 


def get_read_ids(barcode, fastq_file, min_distance):

    barcodes = barcode.split("_")
    read_ids_list = []

    count = 0 
    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if (i - 1) % 4 == 0:
                read_id = line
            elif (i - 2) % 4 == 0:
                read = line
                if (
                    hamming_distance(barcodes[0], read[:8]) <= min_distance and
                    hamming_distance(barcodes[1], read[12:20]) <= min_distance and
                    hamming_distance(barcodes[2], read[24:32]) <= min_distance and
                    hamming_distance(barcodes[3], read[36:44]) <= min_distance
                ):
                    read_id = read_id.split()[0]
                    read_id = read_id.replace("@", "")

                    read_ids_list.append(read_id)

    return read_ids_list


def parse_arguments():
    parser = argparse.ArgumentParser(description="Select fastq read IDs based on the dominant barcodes")
    parser.add_argument("--barcode", required=True, help="Barcode where you want to select representative reads for")
    parser.add_argument("--fastq", required=True, help="Path to the fastq file you want to select the reads from")
    parser.add_argument("--outputdir", required=True, help="Directory where to write the sequence read IDs. Filename will consist of the barcode + txt")
    parser.add_argument("--distance", required=False, type=int, default=1, help="maximum hamming distance per barcode of 8 nt")

    return parser.parse_args()


def main():
    args = parse_arguments()

    barcode = args.barcode
    fastq_file = args.fastq
    output_dir = args.outputdir
    mindist = args.distance

    read_ids = get_read_ids(barcode, fastq_file, mindist)
    
    output_file = f"{output_dir}/{barcode}.readids.txt"
    with open(output_file, "w") as f:
        f.writelines(f"{item}\n" for item in read_ids)

 
if __name__ == "__main__":
    main()


# ./demux_barcodes_getreads.py --outputdir /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/results/demux --fastq /user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/data/fastq/scDNA_AT_01_2.fq.gz --barcode GGTCTCAT_AACTGCTC_CACCATCT_GTGACTCT