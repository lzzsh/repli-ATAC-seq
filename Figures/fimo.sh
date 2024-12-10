#!/bin/bash
#SBATCH -c 30
#SBATCH --mem=128GB
#SBATCH -J Fimo
#SBATCH -o Fimo.%j.out
#SBATCH -e Fimo.%j.err
#SBATCH -p lxd_20241204,intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

source /home/liuxiaodongLab/liaozizhuo/anaconda3/etc/profile.d/conda.sh
conda activate py37
fimo --thresh 1e-3 -oc meme JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta 