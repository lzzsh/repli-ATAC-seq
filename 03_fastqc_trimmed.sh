#!/bin/bash
#SBATCH -c 6
#SBATCH --mem=16GB
#SBATCH -J fastq
#SBATCH -o fastq.%j.out
#SBATCH -e fastq.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load /soft/modules/modulefiles/bioinfo/fastqc/0.11.9

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/
mkdir fastqc_trimmed_rawdata

for file in $(ls ./trimmomatic_rawdata/*P.fq.gz);do
	fastqc ${file} -o ./fastqc_trimmed_rawdata
	sleep 1
done
