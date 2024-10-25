#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J bwa
#SBATCH -o bwa.%j.out
#SBATCH -e bwa.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load samtools/1.13
module load bwa/0.7.17

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/trimmomatic_rawdata
mkdir ../bwa_all_rawdata

out="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata"
for file in $(ls *1P.fq.gz);do
	prefix=${file%%_*}
	bwa mem rice_all_genomes_v7_bwa_index -t 10 -M -k 32 -R '@RG\tID:${prefix}\tSM:${prefix}\tLB:${prefix}\tPL:illumina' ${prefix}_1P.fq.gz ${prefix}_2P.fq.gz | samtools view -ubhS > ${out}/${prefix}.bam
	sleep 1
done

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata

for f in $(ls *bam);do
	prefix=${f%%\.*}
	samtools flagstat $f > ${prefix}_flagstat.out
done