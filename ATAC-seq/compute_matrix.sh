## WT
computeMatrix reference-point  -b 2000 -a 2000 -R /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/all_ES.bed \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/all_MS.bed \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/all_LS.bed -S \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/bws/LZZ-18-1.bw \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/bws/LZZ-18-4.bw \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/bws/LZZ-18-2.bw \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/bws/LZZ-18-3.bw \
--skipZeros --missingDataAsZero -o matrix_all_peaks_WT.gz --outFileSortedRegions regions_all_peaks_WT.bed --outFileNameMatrix matrix_all_peaks_WT.tsv.gz

plotProfile -m matrix_all_peaks_WT.gz  -out TSS_profile_all_peaks_WT.pdf  --refPointLabel peak --regionsLabel ES MS LS --outFileSortedRegions location_all_peaks_WT.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_all_peaks_WT.gz  -out TSS_heatmap_all_peaks_WT.pdf  --refPointLabel peak --regionsLabel ES MS LS --outFileSortedRegions location_all_peaks_WT.bed --plotFileFormat pdf --perGroup \
--samplesLabel signal_ES-1 signal_ES-2 signal_MS signal_LS --heatmapHeight 20 --heatmapWidth 5 --colorMap OrRd 

## TCX2-1 CR
computeMatrix reference-point  -b 2000 -a 2000 -R /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/all_ES.bed \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/all_MS.bed \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/all_LS.bed -S \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/bws/LZZ-18-5.bw \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/bws/LZZ-18-8.bw \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/bws/LZZ-18-6.bw \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/bwa_all_rawdata/bws/LZZ-18-7.bw \
--skipZeros --missingDataAsZero -o matrix_all_peaks_TCX2_1.gz --outFileSortedRegions regions_all_peaks_TCX2_1.bed --outFileNameMatrix matrix_all_peaks_TCX2_1.tsv.gz

plotProfile -m matrix_all_peaks_TCX2_1.gz  -out TSS_profile_all_peaks_TCX2_1.pdf  --refPointLabel peak --regionsLabel ES MS LS --outFileSortedRegions location_all_peaks_TCX2_1.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_all_peaks_TCX2_1.gz  -out TSS_heatmap_all_peaks_TCX2_1.pdf  --refPointLabel peak --regionsLabel ES MS LS --outFileSortedRegions location_all_peaks_TCX2_1.bed --plotFileFormat pdf --perGroup \
--samplesLabel signal_ES-1 signal_ES-2 signal_MS signal_LS --heatmapHeight 20 --heatmapWidth 5 --colorMap OrRd 
