#!/usr/bin/env python3

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=0:50:00
#SBATCH --job-name=count
#SBATCH -A ap_itg_mpu


import gzip
from collections import defaultdict
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process FASTQ files to count barcodes.")
    parser.add_argument("--fastq", required=True, help="Path to the input FASTQ file")
    parser.add_argument("--output", help="Path to the output barcode count file. If not provided, inferred from FASTQ file.")
    parser.add_argument("--limit", type=int, default=None, help="Optional limit to the number of reads to process.")
    parser.add_argument("--mincount", type=int, default=10000, help="Optional number to limit the number of cells to only those cells that have more than this number of barcodes")
    parser.add_argument("--whitelist", default="/user/antwerpen/205/vsc20587/scratch/leishmania_scDNA_atrandi/data/atrandi/atrandi_whitelist_barcodes.csv", help="Put a whitelist file of Atrandi here, or use the default")

    return parser.parse_args()

def process_fastq(fastq_file, limit=None):
    count_barcode = defaultdict(int)
    count = 0

    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f, start=1):
            if (i - 2) % 4 == 0:  # Process the sequence line
                count += 1
                read = line.strip()

                barcode = f"{read[:8]}_{read[12:20]}_{read[24:32]}_{read[36:44]}"
                if 'N' not in barcode:
                    count_barcode[barcode] += 1

                if limit and count >= limit:
                    break

    return count_barcode


def read_whitelist(whitelist_file):
    ## create empty dictionary
    lists_by_pos = defaultdict(list)

    # Open the file and process the data
    with open(whitelist_file, 'r') as file:
        # Skip the header line
        next(file)
        
        # Process each line in the file
        for line in file:
            columns = line.strip().split()
            pos, seq = columns[0], columns[2]
            lists_by_pos[pos].append(seq)

    return lists_by_pos


def check_barcode(barcode, whitelist_barcode):
    barcodes = barcode.split("_")
    for i, segment in enumerate(barcodes, start=1):
        position = str(i)  # Convert position to string for dictionary lookup
        if segment not in whitelist_barcode[str(4-i+1)]:
            print(f"barcode {segment} is not part of whitelist barcodes at position {i}")
            print(whitelist_barcode[str(4-i+1)])
            return False  # Segment not found in the corresponding whitelist
    return True


def write_barcode_counts(barcode_counts, output_file, mincount, whitelist_file):
    sorted_barcodes = sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True)

    print(whitelist_file)
    whitelist_barcode = read_whitelist(whitelist_file)

    with open(output_file, 'w') as out_fh:
        for barcode, count in sorted_barcodes:
            
            ## break out to only keep the most dominant barcodes
            if count < mincount:
                break

            ## check whether the barcodes are part of the whitelist
            if not check_barcode(barcode, whitelist_barcode):
                print(f"barcode {barcode} is not compliant with whitelist barcodes")
                continue

            line_out = f"{barcode}\t{count}\n"
            print(line_out)
            out_fh.write(line_out)

def main():
    args = parse_arguments()

    fastq_file = args.fastq
    output_file = args.output or fastq_file.replace("_2.fq.gz", "_barcode_counts.txt")
    
    print(f"Processing FASTQ file: {fastq_file}")
    if args.limit:
        print(f"Processing will stop after {args.limit} reads.")

    mincount = args.mincount
    if mincount:
        print(f"Only cells with more than {mincount} reads will be processed.")

    whitelist_file = args.whitelist

    barcode_counts = process_fastq(fastq_file, limit=args.limit)

    write_barcode_counts(barcode_counts, output_file, mincount, whitelist_file)
    print(f"barcode counts written to: {output_file}")

if __name__ == "__main__":
    main()
