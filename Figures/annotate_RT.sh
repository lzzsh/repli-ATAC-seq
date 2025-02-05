#!/bin/bash

# IDR_peaks
bedtools intersect -a ZH11_RT.gff3 -b annotate.bed -wa -wb > annotate_RT.bed
cut -f 1,2,3,4,8 annotate_RT.bed > annotate_RT_temp.bed
sort -k2,2n annotate_RT_temp.bed | uniq > annotate_RT.bed

# ZH11-2_peaks
bedtools intersect -a ZH11_2_RT_control.gff3 -b annotate.bed -wa -wb > annotate_RT_2.bed
cut -f 1,2,3,4,8 annotate_RT_2.bed > annotate_RT_temp.bed 
mv annotate_RT_temp.bed annotate_RT_2.bed