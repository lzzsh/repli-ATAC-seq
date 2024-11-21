#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=16GB
#SBATCH -J trim
#SBATCH -o trim.%j.out
#SBATCH -e trim.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load trimmomatic/0.39

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq

# Create output directory
mkdir -p trimmomatic_rawdata

# Function to rename files
rename_files() {
    for file in ./01.RawData/*R1.fastq.gz; do
        new_file=$(echo $file | sed 's/R1/1/g' | sed 's/fastq/fq/g')
        mv "$file" "$new_file"
    done

    for file in ./01.RawData/*R2.fastq.gz; do
        new_file=$(echo $file | sed 's/R2/2/g' | sed 's/fastq/fq/g')
        mv "$file" "$new_file"
    done
}

# Function to run Trimmomatic
run_trimmomatic() {
    for file in $(ls ./01.RawData/*1.fq.gz); do
        file=$(basename $file)
        prefix=${file%_*}
        trimmomatic PE -threads 10 -phred33 -trimlog ./trimmomatic_rawdata/${prefix}trim.log \
            -basein ./01.RawData/${prefix}_1.fq.gz -baseout ./trimmomatic_rawdata/${prefix}.fq.gz \
            ILLUMINACLIP:/storage/publicdata/trimmomatic_adapters/NexteraPE-PE.fa:2:30:10:8:TRUE \
            SLIDINGWINDOW:4:15 MINLEN:30
        sleep 1
    done
}

# Execute renaming and Trimmomatic
rename_files
run_trimmomatic