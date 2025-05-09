bedtools makewindows -g /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/genome.chrom.sizes.sorted -w 1000 -s 500 > rice_1kb_step500_windows.bed


bedtools multicov -bams \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-2_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-3_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-4_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-5_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-6_sort_rmdup_rmor_q30.bam \
-bed rice_1kb_step500_windows.bed  > raw_counts.txt