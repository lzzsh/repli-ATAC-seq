# motif segmentation overlap
awk '$4=="MA1379.1" {print $1"\t"$2"\t"$3}' fimo.txt > SOL1_MA1379.1.bed
awk '$4=="MA1162.1" {print $1"\t"$2"\t"$3}' fimo.txt > TCX2_MA1162.1.bed

bedtools slop -i SOL1_MA1379.1.bed -g /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/genome.chrom.sizes.sorted -b 500 > SOL1_peaks_1kb.bed
bedtools slop -i TCX2_MA1162.1.bed -g /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/genome.chrom.sizes.sorted -b 500 > TCX2_peaks_1kb.bed

bedtools intersect -a segmentation_all.bed -b SOL1_peaks_1kb.bed -wa |sort|uniq > SOL1_overlap_segmentation.bed
bedtools intersect -a segmentation_all.bed -b TCX2_peaks_1kb.bed -wa |sort|uniq > TCX2_overlap_segmentation.bed

# motif genome overlap
bedtools intersect -a segmentation_all.bed -b SOL1_MA1379.1.bed -wa |sort|uniq > SOL1_overlap_genome.bed
bedtools intersect -a segmentation_all.bed -b TCX2_MA1162.1.bed -wa |sort|uniq > TCX2_overlap_genome.bed


# open chromatin
# motif genome overlap
bedtools intersect -a /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT_all_org.gff3 -b SOL1_MA1379.1.bed -wa |sort|uniq > SOL1_overlap_peak.bed
bedtools intersect -a /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT_all_org.gff3 -b TCX2_MA1162.1.bed -wa |sort|uniq > TCX2_overlap_peak.bed