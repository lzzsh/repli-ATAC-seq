#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J featureCounts
#SBATCH -o featureCounts.%j.out
#SBATCH -e featureCounts.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load subread/2.0.3   

BAM_DIR="/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/star_all_rawdata"
GTF="/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3"
OUT_DIR="${BAM_DIR}/featureCounts_out"
mkdir -p ${OUT_DIR}

BAM_FILES=$(ls ${BAM_DIR}/*_Aligned.sortedByCoord.out.bam)

featureCounts -T 8 -p -B -C \
  -a ${GTF} \
  -o ${OUT_DIR}/gene_counts.txt \
  -g gene_id \
  -t exon \
  ${BAM_FILES}


echo "featureCounts finished successfully!"
echo "Output saved to: ${OUT_DIR}/gene_counts.txt"