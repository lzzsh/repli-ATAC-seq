#!/bin/bash

# NIP and ZH11-1
# select the peaks with IDR < 0.05
awk '{if($5 >= 540) print $0}' ES-idr > ES-idr_0.05
awk '{if($5 >= 540) print $0}' MS-idr > MS-idr_0.05
awk '{if($5 >= 540) print $0}' LS-idr > LS-idr_0.05

# sort and merge peaks
cat ES-idr_0.05 MS-idr_0.05 LS-idr_0.05 > EdU-idr_0.05
sort -k1,1 -k2,2n EdU-idr_0.05 > EdU-idr_0.05_sort
bedtools merge -i EdU-idr_0.05_sort > EdU-idr_0.05_sort_merge

# calculate reads number in each open chromatin region
location="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData"
bedtools multicov -bams \
${location}/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam \
${location}/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam \
${location}/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam \
-bed EdU-idr_0.05_sort_merge > EdU-idr_reads_ZH11.txt

# select the peaks with IDR < 0.1
awk '{if($5 >= 415) print $0}' ES-idr > ES-idr_0.1
awk '{if($5 >= 415) print $0}' MS-idr > MS-idr_0.1
awk '{if($5 >= 415) print $0}' LS-idr > LS-idr_0.1

# sort and merge peaks
cat ES-idr_0.1 MS-idr_0.1 LS-idr_0.1 > EdU-idr_0.1
sort -k1,1 -k2,2n EdU-idr_0.1 > EdU-idr_0.1_sort
bedtools merge -i EdU-idr_0.1_sort > EdU-idr_0.1_control_sort_merge

# calculate reads number in each open chromatin region
location="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData"
bedtools multicov -bams \
${location}/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam \
${location}/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam \
${location}/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam \
-bed EdU-idr_0.1_control_sort_merge > EdU-idr_reads_ZH11_control.txt


# ZH11-2 (no replicate)
cat /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_p0.01/ZH11_2_LS_out/ZH11_2_LS_peaks_sort.narrowPeak \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_p0.01/ZH11_2_ES_out/ZH11_2_ES_peaks_sort.narrowPeak \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_p0.01/ZH11_2_MS_out/ZH11_2_MS_peaks_sort.narrowPeak > ZH11_0.1
sort -k1,1 -k2,2n ZH11_0.1 > ZH11_0.1_sort
bedtools merge -i ZH11_0.1_sort > ZH11_0.1_control_sort_merge
bedtools multicov -bams /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam \
	-bed ZH11_0.1_control_sort_merge > ZH11_reads_ZH11-2_control.txt

# wox11
cat /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/macs2/macs2_p0.01/xw11_CR_LS_out/xw11_CR_LS_peaks_sort.narrowPeak \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/macs2/macs2_p0.01/xw11_CR_ES_out/xw11_CR_ES_peaks_sort.narrowPeak \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/macs2/macs2_p0.01/xw11_CR_MS_out/xw11_CR_MS_peaks_sort.narrowPeak > xw11_CR_0.1
sort -k1,1 -k2,2n xw11_CR_0.1 > xw11_CR_0.1_sort
bedtools merge -i xw11_CR_0.1_sort > xw11_CR_0.1_sort_merge
bedtools multicov -bams /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/LZZ-2-g_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-h_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-i_S99_sort_rmdup_rmor_q30.bam \
	-bed xw11_CR_0.1_sort_merge > xw11_CR_reads.txt

# x39
cat /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/macs2/macs2_p0.01/x39_CR_LS_out/x39_CR_LS_peaks_sort.narrowPeak \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/macs2/macs2_p0.01/x39_CR_ES_out/x39_CR_ES_peaks_sort.narrowPeak \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/macs2/macs2_p0.01/x39_CR_MS_out/x39_CR_MS_peaks_sort.narrowPeak > x39_CR_0.1
sort -k1,1 -k2,2n x39_CR_0.1 > x39_CR_0.1_sort
bedtools merge -i x39_CR_0.1_sort > x39_CR_0.1_sort_merge
bedtools multicov -bams /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/LZZ-3-j_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-k_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-l_S99_sort_rmdup_rmor_q30.bam \
	-bed x39_CR_0.1_sort_merge > x39_CR_reads.txt

