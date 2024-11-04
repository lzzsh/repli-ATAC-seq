#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J rmor
#SBATCH -o rmor.%j.out
#SBATCH -e rmor.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load samtools/1.13

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata

for file in $(ls *_sort_rmdup.bam);do
	prefix=${file%\.*}
	samtools view -h ${file} | grep -v 'chr13\|chr14\|chrSy\|chrUn' | samtools view -bh > ${prefix}_rmor.bam
	sleep 1
done