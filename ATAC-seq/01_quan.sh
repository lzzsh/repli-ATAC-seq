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
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-19_LZZ-19-1_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-19_LZZ-19-4_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-19_LZZ-19-2_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-19_LZZ-19-3_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-20_LZZ-20-9_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-20_LZZ-20-12_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-20_LZZ-20-10_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-20_LZZ-20-11_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-20_LZZ-20-17_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-20_LZZ-20-20_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-20_LZZ-20-18_sort_rmdup_rmor_q30.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/bwa_all_rawdata/LZZ-20_LZZ-20-19_sort_rmdup_rmor_q30.bam \
-bed /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/peaks_quan_tpm_filtered_pos_org.txt  > ATAC_peaks_quan_peak.txt