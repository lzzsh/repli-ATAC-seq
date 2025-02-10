awk '{print $1"\t"$2"\t"$3}' peaks_quan_tpm_filtered.txt > peaks_quan_tpm_filtered_pos.txt

# select reads only from open chromatin region
bedtools multicov -bams \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-2_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-3_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-d-sup_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-1_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-3_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-6_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-4_sort_rmdup_rmor_q30.bam \
	/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-5_sort_rmdup_rmor_q30.bam \
-bed peaks_quan_tpm_filtered_pos.txt  > peaks_re_quan.txt

# select reads only from open chromatin region
files=(
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-2_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-3_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-p_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-j_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/LZZ-2-g_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-d-sup_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-1_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-3_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-4_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-5_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-6_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-q_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-k_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-n_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-h_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-r_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-l_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-o_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-i_S99_sort_rmdup_rmor_q30.bam"
	)

for file in ${files[@]}; do
	name=$(basename ${file})
	name=${name%%.*}
	bedtools intersect -abam ${file} -b peaks_quan_tpm_filtered_pos.txt -wa > ${name}_filtered.bam
    samtools index ${name}_filtered.bam
done

#mkdir bws
for file in $(find ./ -name '*_filtered.bam');do
	name=$(basename ${file})
	prefix=${name%_*}
	bamCoverage --bam ${file} -o ./bws/${prefix}.bw --normalizeUsing BPM -p 5 --binSize 50
	sleep 1
done