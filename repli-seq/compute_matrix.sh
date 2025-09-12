computeMatrix reference-point  -b 2000 -a 2000 -R /storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/all_ES.bed \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/all_MS.bed \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/all_LS.bed -S \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/bws/WT-G1.bw \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/bws/WT-ES.bw \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/bws/LZZ-14-2.bw \
	/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/bwa_all_rawdata/bws/LZZ-16-3.bw \
--skipZeros --missingDataAsZero -o matrix_all_peaks_WT.gz --outFileSortedRegions regions_all_peaks_WT.bed --outFileNameMatrix matrix_all_peaks_WT.tsv.gz

plotProfile -m matrix_all_peaks_WT.gz  -out TSS_profile_all_peaks_WT.pdf  --refPointLabel peak --regionsLabel ES MS LS --outFileSortedRegions location_all_peaks_WT.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_all_peaks_WT.gz  -out TSS_heatmap_all_peaks_WT.pdf  --refPointLabel peak --regionsLabel ES MS LS --outFileSortedRegions location_all_peaks_WT.bed --plotFileFormat pdf --perGroup \
--samplesLabel signal_G1 signal_ES signal_MS signal_LS --heatmapHeight 20 --heatmapWidth 5 --colorMap OrRd 