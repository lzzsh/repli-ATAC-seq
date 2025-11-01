#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J star_align
#SBATCH -o star.%j.out
#SBATCH -e star.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load star/2.7.11b
module load samtools/1.13

RAW_DIR="/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/01.RawData"
REF_FASTA="/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta"
REF_GTF="/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3"
STAR_INDEX="/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/star_index"
OUT_DIR="/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata"

mkdir -p ${OUT_DIR}

if [ ! -d "${STAR_INDEX}" ] || [ -z "$(ls -A ${STAR_INDEX})" ]; then
    echo "STAR genome index not found — building now..."
    mkdir -p ${STAR_INDEX}

    STAR \
        --runThreadN 8 \
        --runMode genomeGenerate \
        --genomeDir ${STAR_INDEX} \
        --genomeFastaFiles ${REF_FASTA} \
        --sjdbGTFfile ${REF_GTF} \
        --sjdbOverhang 99 \
    --genomeSAindexNbases 13

    echo "Genome index built successfully!"
else
    echo "STAR genome index already exists, skipping build."
fi

cd ${RAW_DIR}
for file in $(ls *_R1.fq.gz); do
    prefix=${file%%_R1.fq.gz}
    echo "Processing sample: ${prefix}"

    STAR \
        --runThreadN 8 \
        --genomeDir ${STAR_INDEX} \
        --readFilesIn ${prefix}_R1.fq.gz ${prefix}_R2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${OUT_DIR}/${prefix}_ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMattrRGline ID:${prefix} SM:${prefix} LB:${prefix} PL:illumina

    samtools flagstat ${OUT_DIR}/${prefix}_Aligned.sortedByCoord.out.bam > ${OUT_DIR}/${prefix}_flagstat.out
    sleep 1
done