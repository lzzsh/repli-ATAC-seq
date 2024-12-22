#!/bin/bash
awk '$3 == "gene"' all_DIY.gff3|awk 'BEGIN{FS="\t|=|;";OFS="\t"}{print $1,$4-1,$4,$3,$6,$7}'>TSS.bed

# IDR_peaks
bedtools intersect -a TSS.bed -b ZH11_RT.gff3 -wa -wb > TSS_RT.bed
cut -f 1,2,3,4,5,6,10 TSS_RT.bed > TSS_RT_temp.bed
mv TSS_RT_temp.bed TSS_RT.bed

# Control_IDR_peaks
bedtools intersect -a TSS.bed -b ZH11_RT_control.gff3 -wa -wb > TSS_RT_control.bed
cut -f 1,2,3,4,5,6,10 TSS_RT_control.bed > TSS_RT_control_temp.bed
mv TSS_RT_control_temp.bed TSS_RT_control.bed

# ZH11-2_peaks
bedtools intersect -a TSS.bed -b ZH11_2_RT_control.gff3 -wa -wb > TSS_RT_2_control.bed
cut -f 1,2,3,4,5,6,10 TSS_RT_2_control.bed > TSS_RT_2_control_temp.bed
mv TSS_RT_2_control_temp.bed TSS_RT_2_control.bed

# x39 cuttag
bedtools intersect -a TSS_RT_2_control.bed -b x39_cut_peaks_sort.narrowPeak -wa -wb > x39_RT.bed
cut -f 1,2,3,4,5,6,7 x39_RT.bed > x39_temp.bed
mv x39_temp.bed x39_RT.bed

# x49 cuttag
bedtools intersect -a TSS_RT_2_control.bed -b x49_cut_peaks_sort.narrowPeak -wa -wb > x49_RT.bed
cut -f 1,2,3,4,5,6,7 x49_RT.bed > x49_temp.bed
mv x49_temp.bed x49_RT.bed

# xw11 cuttag
bedtools intersect -a TSS_RT_2_control.bed -b xw11_cut_peaks_sort.narrowPeak -wa -wb > xw11_RT.bed
cut -f 1,2,3,4,5,6,7 xw11_RT.bed > xw11_temp.bed
mv xw11_temp.bed xw11_RT.bed

# top3000

# x39 cuttag
bedtools intersect -a /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/TSS_RT_2_control.bed -b ../x39_cut_out/x39_cut_peaks_top3000.narrowPeak -wa -wb > ../x39_cut_out/x39_RT_3000.bed
cut -f 1,2,3,4,5,6,7 ../x39_cut_out/x39_RT_3000.bed > ../x39_cut_out/x39_temp.bed
mv ../x39_cut_out/x39_temp.bed ../x39_cut_out/x39_RT_3000.bed

# x49 cuttag
bedtools intersect -a /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/TSS_RT_2_control.bed -b ../x49_cut_out/x49_cut_peaks_top3000.narrowPeak -wa -wb > ../x49_cut_out/x49_RT_3000.bed
cut -f 1,2,3,4,5,6,7 ../x49_cut_out/x49_RT_3000.bed > ../x49_cut_out/x49_temp.bed
mv ../x49_cut_out/x49_temp.bed ../x49_cut_out/x49_RT_3000.bed

# xw11 cuttag
bedtools intersect -a /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/TSS_RT_2_control.bed -b ../xw11_cut_out/xw11_cut_peaks_top3000.narrowPeak -wa -wb > ../xw11_cut_out/xw11_RT_3000.bed
cut -f 1,2,3,4,5,6,7 ../xw11_cut_out/xw11_RT_3000.bed > ../xw11_cut_out/xw11_temp.bed
mv ../xw11_cut_out/xw11_temp.bed ../xw11_cut_out/xw11_RT_3000.bed