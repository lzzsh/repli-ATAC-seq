#!/bin/bash
bedtools intersect -a ZH11_RT.gff3 -b annotate.bed -wa -wb > annotate_RT.bed
cut -f 1,2,3,4,8 annotate_RT.bed > annotate_RT_temp.bed 
mv annotate_RT_temp.bed annotate_RT.bed