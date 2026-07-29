#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=WGS_GCOX
#SBATCH --mail-type=BEGIN,END,FAIL

#This is the script used to perform the analysis on the bulk whole genome sequencing (WGS) data of both HU3 and BPK081 strains.
#This script is meant to be ran in the VSC HPC at the root folder of the `2025_scDNA_paper` repository.

#prepare the environment
set -euo pipefail
IFS=$'\n\t'


#-----code start---------
##set the directories
REPO_ROOT="$VSC_SCRATCH/projects/2025_scDNA_paper"
outputs_dir="$REPO_ROOT/inputs/WGS_data" #must be absolute path

##get the global variables and functions
source "$REPO_ROOT/set_env.sh"
##create the directory to store the data
rm -rf "$outputs_dir"
mkdir -p "$outputs_dir/raw_fastq_files"
cd "$outputs_dir/raw_fastq_files"
##copy the files
echo copying fastq files to inputs/WGS...
find "/data/antwerpen/grp/aitg/jcdujardin/LeishMAP_WGS_20260703_gabriel_GC-OX-MFA_samples" -type f -name "MFA_4*.fq.gz" -exec cp {} . \;
find "/data/antwerpen/grp/aitg/jcdujardin/LeishMAP_WGS_20260703_gabriel_GC-OX-MFA_samples" -type f -name "MFA_6*.fq.gz" -exec cp {} . \;

#rename files to match the requirements of the Illumina_tools and WGS_tools modules
##requirements are:
## 1) file extension must be `fastq.gz`, not `fq.gz`
## 2) read1 file must end with `R1.fastq.gz`
## 3) read2 file must end with `R2.fastq.gz`
## 4) There must be 1 file pair per sample (so multiple lanes/runs must be merged)!
## 5) File names should be <sample name>_R1.fastq.gz (for read1) or <sample name>_R2.fastq.gz (for read2) only. 
## 6) Sample name cannot contain underscores(_).
echo renaming files...
rename -v '_1.fq.gz' '_R1.fastq.gz' *.fq.gz
rename -v '_2.fq.gz' '_R2.fastq.gz' *.fq.gz
rename -v 'MFA_4' 'BPK081' *.fastq.gz
rename -v 'MFA_6' 'HU3' *.fastq.gz

#merge the files
echo "merging files from multiple lanes into a single file"
    for sample in $(ls *.fastq.gz 2>/dev/null | cut -d'_' -f1 | sort -u); do
        r1_files=$(ls ${sample}_*_R1.fastq.gz | sort -V)
        r2_files=$(ls ${sample}_*_R2.fastq.gz | sort -V)

        echo "Sample: $sample"
        echo "  R1 -> ${sample}_R1.fastq.gz"
        printf '    %s\n' $r1_files
        echo "  R2 -> ${sample}_R2.fastq.gz"
        printf '    %s\n' $r2_files
        echo

        cat $r1_files > ${sample}_R1.fastq.gz &&
        cat $r2_files > ${sample}_R2.fastq.gz &&
        rm $r1_files &&
        rm $r2_files
    done

#----- Run the WGS pipeline ---------
#now that file names are clean, we write a tsv file listing them to be able to submit them as slurm arrays.
cd "$outputs_dir"
file_list -i raw_fastq_files -o temp.tsv

#now we use the tsv file to send fastp as a slurm array
echo "cleaning reads with fastp"
nfiles=$( cat temp.tsv | wc -l )
run_fastp \
    -l "temp.tsv" \
    -o "cleaned_fastq_files/"

#now create another tsv file specifying the sample name, read1 and read2 paths for the cleaned reads
rm temp.tsv
file_list -i "cleaned_fastq_files/" -o temp.tsv

#use it to run BWA-MEM 
echo "mapping reads with bwa-mem"
nfiles=$( cat temp.tsv | wc -l )
run_bwa \
    -l temp.tsv \
    -r $WGS_REF \
    -o "bam_files/"

#now use the generated bam files to count binned reads
# Should be:
echo "calculating binned counts..."
find "bam_files/" -maxdepth 1 -name "*.bam" > temp.tsv
nfiles=$( cat temp.tsv | wc -l )
binned_counts \
    -l temp.tsv \
    -r "$WGS_REF" \
    -o "count_files/" \
    -s $WGS_BIN_SIZE

#finally we estimate genome mappability and bind it to the output count files
echo "estimating bin mappability..."
rm -f *mappability.tsv
ref_mappability \
    --ref "$WGS_REF" \
    --output-dir . \
    --bin-size $WGS_BIN_SIZE

map_file=$( echo *mappability.tsv )

#then we merge the files
tail -n +2 "$map_file" | cut -f9 > mappability_column.temp

for bed_file in "${outputs_dir}/count_files/"*.bed; do
    paste "$bed_file" mappability_column.temp > "${bed_file}.temp"
    mv "${bed_file}.temp" "$bed_file"
done

#remove temporary files
rm *.temp

#now perform aneuploidy analysis
echo "estimating aneuploidy"
get_aneuploidy \
    --inputs_dir="count_files/" \
    --outputs_dir="aneuploidy_analysis/" \
    --ploidy=2

echo "[DONE] Files stored in ${outputs_dir}/"
