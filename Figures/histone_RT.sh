#!/bin/bash

# IDR_peaks
cut -f 1,2,3,9,10 Histone.bed > Histone_temp.bed
mv Histone_temp.bed Histone.bed 
awk 'BEGIN {FS=OFS="\t"} {if ($1 ~ /^chr[1-9]$/) $1="chr0" substr($1, 4); print}' Histone.bed > Histone_temp.bed
mv Histone_temp.bed Histone.bed 
bedtools intersect -a ZH11_RT.gff3 -b Histone.bed -wa -wb > Histone_RT.bed

# ZH11-2_peaks
bedtools intersect -a ZH11_2_RT_control.gff3 -b Histone.bed -wa -wb > Histone_RT_2.bed