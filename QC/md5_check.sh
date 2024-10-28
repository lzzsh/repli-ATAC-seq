#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8GB
#SBATCH -J md5
#SBATCH -o md5.%j.out
#SBATCH -e md5.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData

for file in *R1.fastq.gz.md5; do
    prefix=${file%_*}
    echo "Checking ${prefix}_R1.fastq.gz.md5" >> checksum_results.txt
    md5sum -c "${prefix}_R1.fastq.gz.md5" >> checksum_results.txt
    echo "-----------------------------" >> checksum_results.txt
    echo "Checking ${prefix}_R2.fastq.gz.md5" >> checksum_results.txt
    md5sum -c "${prefix}_R2.fastq.gz.md5" >> checksum_results.txt
    echo "-----------------------------" >> checksum_results.txt
done