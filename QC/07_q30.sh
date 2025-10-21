#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J q30
#SBATCH -o q30.%j.out
#SBATCH -e q30.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load samtools/1.13

cd /storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata

for file in $(ls *_rmor.bam);do
	prefix=${file%\.*}
	samtools view -bq 30 ${file} | samtools view -bF 4 -F 256 > ${prefix}_q30.bam 
	sleep 1
done