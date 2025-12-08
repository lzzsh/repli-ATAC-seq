# IDR_peaks
bedtools intersect -a /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/TSS.bed -b /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation/ZH11_RT_bin.bed -wa -wb > TSS_RT_bin.bed
cut -f 1,2,3,4,5,6,10,11,12,13,14 TSS_RT_bin.bed > TSS_RT_temp.bed
mv TSS_RT_temp.bed TSS_RT_bin.bed
