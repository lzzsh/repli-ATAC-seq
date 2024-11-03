#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J tss
#SBATCH -o tss.%j.out
#SBATCH -e tss.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load samtools/1.13
module load R/4.4.0

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/
mkdir tss_plot
Rscript tss_enrichment.R