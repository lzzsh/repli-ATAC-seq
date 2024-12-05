#!/bin/bash

# IDR_peaks
# Sort the files
sort -k1,1 -k2,2n ZH11_RT.gff3 > ZH11_RT_sorted.gff3
sort -k1,1 -k2,2n annotate_RT.bed > annotate_RT_sorted.bed

# Use bedtools subtract to find unannotated regions
bedtools intersect -a ZH11_RT_sorted.gff3 -b annotate_RT_sorted.bed -v > unannotated_RT.bed

# Delete the intermediate files
rm ZH11_RT_sorted.gff3
rm annotate_RT_sorted.bed


# ZH11-2_peaks
# Sort the files
sort -k1,1 -k2,2n ZH11_2_RT_control.gff3 > ZH11_2_RT_sorted.gff3
sort -k1,1 -k2,2n annotate_RT_2.bed > annotate_RT_2_sorted.bed

# Use bedtools subtract to find unannotated regions
bedtools intersect -a ZH11_2_RT_sorted.gff3 -b annotate_RT_2_sorted.bed -v > unannotated_RT_2.bed

# Delete the intermediate files
rm ZH11_2_RT_sorted.gff3
rm annotate_RT_2_sorted.bed
