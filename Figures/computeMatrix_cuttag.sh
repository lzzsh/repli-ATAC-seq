#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J heatmap
#SBATCH -o heatmap.%j.out
#SBATCH -e heatmap.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn


# ZH11-2_peaks_ES
computeMatrix reference-point  -b 1000 -a 1000 -R /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x39_cut_out/x39_cut_summits.bed \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x49_cut_out/x49_cut_summits.bed \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_cut_summits.bed -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-1-d.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-p.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/bws/LZZ-2-g.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-j.bw  \
--skipZeros -o matrix_peaks_control_CR_2_ES_cut.gz --outFileSortedRegions regions_peaks_control_CR_2_ES_cut.bed

plotProfile -m matrix_peaks_control_CR_2_ES_cut.gz  -out TSS_profile_control_CR_2_ES_cut.pdf  --refPointLabel TF --regionsLabel x39 x49 xw11 --outFileSortedRegions location_control_CR_2_ES_cut.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_control_CR_2_ES_cut.gz  -out TSS_heatmap_control_CR_2_ES_cut.pdf  --refPointLabel TF --regionsLabel x39 x49 xw11 --outFileSortedRegions location_control_CR_2_ES_cut.bed --plotFileFormat pdf --perGroup \
--samplesLabel ES-ZH11-2 ES-wox11-OE ES-wox11-CR ES-x39-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd


# ZH11-2_peaks_MS
computeMatrix reference-point  -b 1000 -a 1000 -R /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x39_cut_out/x39_cut_summits.bed \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x49_cut_out/x49_cut_summits.bed \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_cut_summits.bed -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/bws/LZZ-2-e.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-q.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/bws/LZZ-2-h.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-k.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-n.bw \
--skipZeros -o matrix_peaks_control_CR_2_MS_cut.gz --outFileSortedRegions regions_peaks_control_CR_2_MS_cut.bed

plotProfile -m matrix_peaks_control_CR_2_MS_cut.gz  -out TSS_profile_control_CR_2_MS_cut.pdf  --refPointLabel TF --regionsLabel x39 x49 xw11 --outFileSortedRegions location_control_CR_2_MS_cut.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_control_CR_2_MS_cut.gz  -out TSS_heatmap_control_CR_2_MS_cut.pdf  --refPointLabel TF --regionsLabel x39 x49 xw11 --outFileSortedRegions location_control_CR_2_MS_cut.bed --plotFileFormat pdf --perGroup \
--samplesLabel MS-ZH11-2 MS-wox11-OE MS-wox11-CR MS-x39-CR MS-x49-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd


# ZH11-2_peaks_LS
computeMatrix reference-point  -b 1000 -a 1000 -R /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x39_cut_out/x39_cut_summits.bed \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x49_cut_out/x49_cut_summits.bed \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_cut_summits.bed -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/bws/LZZ-1-f.bw   \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-r.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/bws/LZZ-2-i.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-l.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-o.bw \
--skipZeros -o matrix_peaks_control_CR_2_LS_cut.gz --outFileSortedRegions regions_peaks_control_CR_2_LS_cut.bed

plotProfile -m matrix_peaks_control_CR_2_LS_cut.gz  -out TSS_profile_control_CR_2_LS_cut.pdf  --refPointLabel TF --regionsLabel x39 x49 xw11 --outFileSortedRegions location_control_CR_2_LS_cut.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_control_CR_2_LS_cut.gz  -out TSS_heatmap_control_CR_2_LS_cut.pdf  --refPointLabel TF --regionsLabel x39 x49 xw11 --outFileSortedRegions location_control_CR_2_LS_cut.bed --plotFileFormat pdf --perGroup \
--samplesLabel LS-ZH11-2 LS-wox11-OE LS-wox11-CR LS-x39-CR LS-x49-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd


# x39_top3000_ES
computeMatrix reference-point  -b 1000 -a 1000 -R /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x39_cut_out/x39_ES_peaks_top3000.txt \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x39_cut_out/x39_MS_peaks_top3000.txt \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x39_cut_out/x39_LS_peaks_top3000.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-1-d.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-p.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/bws/LZZ-2-g.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-j.bw  \
--skipZeros -o matrix_peaks_x39_3000_cut.gz --outFileSortedRegions regions_peaks_x39_3000_cut.bed

plotHeatmap -m matrix_peaks_x39_3000_cut.gz  -out TSS_heatmap_x39_3000_cut.pdf  --refPointLabel TF --regionsLabel E M L --outFileSortedRegions location_control_x39_3000_cut.bed --plotFileFormat pdf --perGroup \
--samplesLabel ES-ZH11-2 ES-wox11-OE ES-wox11-CR ES-x39-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd

# x49_top3000_ES
computeMatrix reference-point  -b 1000 -a 1000 -R /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x49_cut_out/x49_ES_peaks_top3000.txt \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x49_cut_out/x49_MS_peaks_top3000.txt \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x49_cut_out/x49_LS_peaks_top3000.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-1-d.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-p.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/bws/LZZ-2-g.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-j.bw  \
--skipZeros -o matrix_peaks_x49_3000_cut.gz --outFileSortedRegions regions_peaks_x49_3000_cut.bed

plotHeatmap -m matrix_peaks_x49_3000_cut.gz  -out TSS_heatmap_x49_3000_cut.pdf  --refPointLabel TF --regionsLabel E M L --outFileSortedRegions location_control_x49_3000_cut.bed --plotFileFormat pdf --perGroup \
--samplesLabel ES-ZH11-2 ES-wox11-OE ES-wox11-CR ES-x39-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd

# xw11_top3000_ES
computeMatrix reference-point  -b 1000 -a 1000 -R /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_ES_peaks_top3000.txt \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_MS_peaks_top3000.txt \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_LS_peaks_top3000.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-1-d.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-p.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/bws/LZZ-2-g.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-j.bw  \
--skipZeros -o matrix_peaks_xw11_3000_cut.gz --outFileSortedRegions regions_peaks_xw11_3000_cut.bed

plotHeatmap -m matrix_peaks_xw11_3000_cut.gz  -out TSS_heatmap_xw11_3000_cut.pdf  --refPointLabel TF --regionsLabel E M L --outFileSortedRegions location_control_xw11_3000_cut.bed --plotFileFormat pdf --perGroup \
--samplesLabel ES-ZH11-2 ES-wox11-OE ES-wox11-CR ES-x39-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd