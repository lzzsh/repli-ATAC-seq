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
computeMatrix reference-point  -b 1000 -a 1000 -R ES_peaks.txt EMS_peaks.txt MS_peaks.txt MLS_peaks.txt LS_peaks.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-1.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-2.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-3.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-4.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-5.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-6.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/bws/LZZ-8-1.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/bws/LZZ-8-2.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/bws/LZZ-8-2.bw \
--skipZeros  -o matrix_peaks.gz --outFileSortedRegions regions_peaks.bed

plotProfile -m matrix_peaks.gz  -out peak_profile.pdf  --refPointLabel peak --regionsLabel E EM M ML L --outFileSortedRegions location.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks.gz  -out peak_heatmap.pdf  --refPointLabel peak --regionsLabel E EM M ML L --outFileSortedRegions location.bed --plotFileFormat pdf --perGroup \
--samplesLabel G1 MS_without_EdU G2 ES MS LS ES_control MS_control LS_control --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd

# Control_IDR_peaks
computeMatrix reference-point  -b 1000 -a 1000 -R ES_control_tss.txt EMS_control_tss.txt MS_control_tss.txt MLS_control_tss.txt LS_control_tss.txt -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-1.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-2.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-3.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-4.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-5.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/bws/origin-2_origin-2-6.bw \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/bws/LZZ-8-1.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/bws/LZZ-8-2.bw /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/bws/LZZ-8-2.bw \
--skipZeros -o matrix_peaks_control.gz --outFileSortedRegions regions_peaks_control.bed

plotProfile -m matrix_peaks_control.gz  -out TSS_profile_control.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_control.gz  -out TSS_heatmap_control.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location_control.bed --plotFileFormat pdf --perGroup \
--samplesLabel G1 MS_without_EdU G2 ES MS LS ES_control MS_control LS_control --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd

plotHeatmap -m matrix_peaks.gz  -out TSS_heatmap.pdf  --refPointLabel TSS --regionsLabel E EM M ML L --outFileSortedRegions location.bed --plotFileFormat pdf --perGroup \
--samplesLabel G1 MS_without_EdU G2 ES MS LS ES_control MS_control LS_control --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd