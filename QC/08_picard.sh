#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J picard
#SBATCH -o picard.%j.out
#SBATCH -e picard.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load picard/2.25.1

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata
mkdir ../insertsize_rawdata

out="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/insertsize_rawdata"
for file in $(ls *_q30.bam);do
	prefix=${file%\.*}
	picard CollectInsertSizeMetrics I=${file} O=${out}/${prefix}_insertsize.txt H=${out}/${prefix}_insertsize.pdf M=0.5 VALIDATION_STRINGENCY=SILENT
	sleep 1
done