#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8GB
#SBATCH -J fastq
#SBATCH -o fastq.%j.out
#SBATCH -e fastq.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load /soft/modules/modulefiles/bioinfo/fastqc/0.11.9

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/
mkdir fastqc_rawdata

for file in $(find ./01.RawData -name '*.fastq.gz');do
	fastqc ${file} -o ./fastqc_rawdata
	sleep 1
done
