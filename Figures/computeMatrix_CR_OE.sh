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

# IDR_peaks
computeMatrix reference-point  -b 1000 -a 1000 -R ES_TSS.txt EMS_TSS.txt MS_TSS.txt MLS_TSS.txt LS_TSS.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-1.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-2.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-3.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-4.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-5.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-6.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/bws/LZZ-8-1.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/bws/LZZ-8-2.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/bws/LZZ-8-2.bw \
--skipZeros  -o matrix_peaks.gz --outFileSortedRegions regions_peaks.bed

plotProfile -m matrix_peaks_CR.gz  -out peak_profile_CR.pdf  --refPointLabel peak --regionsLabel E EM M ML L --outFileSortedRegions location_CR.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_CR.gz  -out peak_heatmap_CR.pdf  --refPointLabel peak --regionsLabel E EM M ML L --outFileSortedRegions location_CR.bed --plotFileFormat pdf --perGroup \
--samplesLabel G1 MS_without_EdU G2 ES MS LS ES_control MS_control LS_control --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd


# Control_IDR_peaks
computeMatrix reference-point  -b 1000 -a 1000 -R ES_control_tss.txt EMS_control_tss.txt MS_control_tss.txt MLS_control_tss.txt LS_control_tss.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/bws/origin-origin-6.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-4.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/bws/LZZ-1-d.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-p.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/bws/LZZ-2-g.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-j.bw  \
--skipZeros -o matrix_peaks_control_CR.gz --outFileSortedRegions regions_peaks_control_CR.bed

plotProfile -m matrix_peaks_control_CR.gz  -out TSS_profile_control_CR.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control_CR.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_control_CR.gz  -out TSS_heatmap_control_CR.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control_CR.bed --plotFileFormat pdf --perGroup \
--samplesLabel ES-NIP ES-ZH11-1 ES-ZH11-2 ES-wox11-OE ES-wox11-CR ES-x39-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd


# ZH11-2_peaks_ES
computeMatrix reference-point  -b 1000 -a 1000 -R ES_2_control_tss.txt EMS_2_control_tss.txt MS_2_control_tss.txt MLS_2_control_tss.txt LS_2_control_tss.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-1-d.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-p.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/bws/LZZ-2-g.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/bws/LZZ-3-j.bw  \
--skipZeros -o matrix_peaks_control_CR_2_ES.gz --outFileSortedRegions regions_peaks_control_CR_2_ES.bed

plotProfile -m matrix_peaks_control_CR_2_ES.gz  -out TSS_profile_control_CR_2_ES.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control_CR_2_ES.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_control_CR_2_ES.gz  -out TSS_heatmap_control_CR_2_ES.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control_CR_2_ES.bed --plotFileFormat pdf --perGroup \
--samplesLabel ES-ZH11-2 ES-wox11-OE ES-wox11-CR ES-x39-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd


# ZH11-2_peaks_MS
computeMatrix reference-point  -b 1000 -a 1000 -R ES_2_control_tss.txt EMS_2_control_tss.txt MS_2_control_tss.txt MLS_2_control_tss.txt LS_2_control_tss.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/bws/LZZ-2-e.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-q.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/bws/LZZ-2-h.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-k.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-n.bw \
--skipZeros -o matrix_peaks_control_CR_2_MS.gz --outFileSortedRegions regions_peaks_control_CR_2_MS.bed

plotProfile -m matrix_peaks_control_CR_2_MS.gz  -out TSS_profile_control_CR_2_MS.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control_CR_2_MS.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_control_CR_2_MS.gz  -out TSS_heatmap_control_CR_2_MS.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control_CR_2_MS.bed --plotFileFormat pdf --perGroup \
--samplesLabel MS-ZH11-2 MS-wox11-OE MS-wox11-CR MS-x39-CR MS-x49-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd


# ZH11-2_peaks_LS
computeMatrix reference-point  -b 1000 -a 1000 -R ES_2_control_tss.txt EMS_2_control_tss.txt MS_2_control_tss.txt MLS_2_control_tss.txt LS_2_control_tss.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/bws/LZZ-1-f.bw   \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-r.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/bws/LZZ-2-i.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-l.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/bws/LZZ-4-o.bw \
--skipZeros -o matrix_peaks_control_CR_2_LS.gz --outFileSortedRegions regions_peaks_control_CR_2_LS.bed

plotProfile -m matrix_peaks_control_CR_2_LS.gz  -out TSS_profile_control_CR_2_LS.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control_CR_2_LS.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_control_CR_2_LS.gz  -out TSS_heatmap_control_CR_2_LS.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control_CR_2_LS.bed --plotFileFormat pdf --perGroup \
--samplesLabel LS-ZH11-2 LS-wox11-OE LS-wox11-CR LS-x39-CR LS-x49-CR --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd