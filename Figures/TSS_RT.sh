#!/bin/bash
awk '$3 == "gene"' all_DIY.gff3|awk 'BEGIN{FS="\t|=|;";OFS="\t"}{print $1,$4-1,$4,$3,$6,$7}'>TSS.bed
bedtools intersect -a TSS.bed -b ZH11_RT.gff3 -wa -wb > TSS_RT.bed
cut -f 1,2,3,4,5,6,10 TSS_RT.bed > TSS_RT_temp.bed
mv TSS_RT_temp.bed TSS_RT.bed