#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J multicov
#SBATCH -o multicov.%j.out
#SBATCH -e multicov.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load bedtools/2.30.0 

# open chromatin
bedtools multicov -bams \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-1_L4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-4_L4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-2_L4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-3_L4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-5_L4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-8_L4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-6_L4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-7_L4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-9_L4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/LZZ-18-10_L4_sort_rmdup_rmor_q30.bam \
-bed /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/peaks_quan_tpm_filtered_pos_org.txt  > ATAC_peaks_quan_peak.txt