#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J reads_count
#SBATCH -o reads_count.%j.out
#SBATCH -e reads_count.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load samtools/1.13

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata

# Find all files ending with q30.bam and process each one
find . -type f -name '*q30.bam' | while read bam_file; do
    bam_file = $(basename $bam_file)
    # Use samtools to count the sequences in the current file
    reads=$(samtools view -c "$bam_file")
    echo "The reads number of $bam_file: $reads" >> reads_count.txt
done