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

# segmentation 1kb bin
bedtools makewindows -g /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/genome.chrom.sizes.sorted -w 1000 > rice_1kb_windows.bed

bedtools multicov -bams \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-1-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-4-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-1_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-4_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-2_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-3_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-5-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-9-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-5_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-9_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-6_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-7_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-10-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-14-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-10_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-14_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-11_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-12_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-15-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-18-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-15_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-18_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-16_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-17_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-19-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-19_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-20_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-21_S99_unique_sort.bam \
-bed rice_1kb_windows.bed  > repli_peaks_quan.txt


# segmentation 10kb bin
bedtools makewindows -g /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/genome.chrom.sizes.sorted -w 10000 > rice_10kb_windows.bed

bedtools multicov -bams \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-1-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-4-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-1_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-4_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-2_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-3_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-5-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-9-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-5_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-9_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-6_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-7_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-10-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-14-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-10_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-14_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-11_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-12_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-15-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-18-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-15_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-18_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-16_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-17_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-19-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-19_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-20_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-21_S99_unique_sort.bam \
-bed rice_10kb_windows.bed  > repli_peaks_quan_10kb.txt

# open chromatin
bedtools multicov -bams \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-1-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-4-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-1_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-4_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-2_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-3_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-5-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-9-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-5_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-9_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-6_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-7_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-10-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-14-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-10_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-14_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-11_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-12_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-15-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-18-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-15_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-18_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-16_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-17_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-13-19-sup_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-12-19_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-14-20_S99_unique_sort.bam \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/LZZ-16-21_S99_unique_sort.bam \
-bed /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/peaks_quan_tpm_filtered_pos_org.txt  > repli_peaks_quan_peak.txt