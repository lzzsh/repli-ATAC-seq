#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J bws
#SBATCH -o bws.%j.out
#SBATCH -e bws.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load samtools/1.13

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata
mkdir bws

#make index
for bam in ./*_q30.bam; do
    samtools index "$bam"
done

#mkdir bws
for file in $(find ./ -name '*q30.bam');do
	name=$(basename ${file})
	prefix=${name%%_*}
	bamCoverage --bam ${file} -o ./bws/${prefix}.bw --normalizeUsing BPM -p 5 --binSize 50
	sleep 1
done