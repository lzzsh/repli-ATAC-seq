awk '{print $1, $2-500, $3+500}' OFS="\t" peaks_quan_tpm_filtered_pos_org.txt > peaks_extended.bed
awk '{if ($2 < 0) $2=0; print $0}' OFS="\t" peaks_extended.bed > peaks_extended_fixed.bed

bedtools multicov -bams \
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
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-p_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-q_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-r_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/LZZ-2-g_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-h_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-i_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-j_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-k_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-l_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-2_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-3_S99_sort_rmdup_rmor_q30.bam \
-bed peaks_extended_fixed.bed > peaks_open_region_counts.txt

awk '{center=int(($2+$3)/2); print $1, center-5000, center+5000}' OFS="\t" peaks_quan_tpm_filtered_pos_org.txt > peaks_5kb_window.bed
awk '{if ($2 < 0) $2=0; print $0}' OFS="\t" peaks_5kb_window.bed > peaks_5kb_window_fixed.bed

bedtools multicov -bams \
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
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-p_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-q_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-r_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/LZZ-2-g_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-h_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-i_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-j_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-k_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-l_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-2_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-3_S99_sort_rmdup_rmor_q30.bam \
-bed peaks_5kb_window_fixed.bed > peaks_5kb_window_counts.txt

paste peaks_5kb_window_counts.txt peaks_open_region_counts.txt | awk '{printf "%s\t%s\t%s", $1, $2, $3; for(i=4; i<=NF/2; i++) printf "\t%d", $i-$(i+NF/2); print ""}' > peaks_masked_counts.txt