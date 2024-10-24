#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32GB
#SBATCH -J rmdup
#SBATCH -o rmdup.%j.out
#SBATCH -e rmdup.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load

cd /public/home/whli/data/ATACellstage/bwa_all_rawdata
for file in $(ls *bam);do
	prefix=${file%%\.*}
	samtools sort -n ${prefix}.bam -o ${prefix}_sort.bam 
	samtools fixmate -cm -O bam ${prefix}_sort.bam - | samtools sort - | samtools markdup -r -s - ${prefix}_sort_rmdup.bam 2>${prefix}_sort_rmdup.log
	sleep 1
done