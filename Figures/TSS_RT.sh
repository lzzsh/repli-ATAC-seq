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