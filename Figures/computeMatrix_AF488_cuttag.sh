computeMatrix reference-point  -b 1000 -a 1000 -R /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x39_cut_out/x39_cut_summits.bed \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/x49_cut_out/x49_cut_summits.bed \
	/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_cut_summits.bed -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/bws/LZZ-8-1_S99_sort_rmdup_rmor_q30.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/bws/LZZ-8-2_S99_sort_rmdup_rmor_q30.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/bws/LZZ-8-3_S99_sort_rmdup_rmor_q30.bw  \
--skipZeros -o matrix_peaks_AF488_cut.gz --outFileSortedRegions regions_peaks_AF488_cut.bed

plotProfile -m matrix_peaks_AF488_cut.gz  -out TSS_profile_AF488_cut.pdf  --refPointLabel TF --regionsLabel x39 x49 xw11 --outFileSortedRegions location_AF488_cut.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_AF488_cut.gz  -out TSS_heatmap_AF488_cut.pdf  --refPointLabel TF --regionsLabel x39 x49 xw11 --outFileSortedRegions location_AF488_cut.bed --plotFileFormat pdf --perGroup \
--samplesLabel ES MS LS --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd

computeMatrix reference-point  -b 1000 -a 1000 -R ES_peaks_center.bed \
	EMS_peaks_center.bed \
	MS_peaks_center.bed \
	LS_peaks_center.bed -S \
/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/bws/x39.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/bws/x49.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/bws/xw11.bw  \
/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/bws/kong.bw \
--skipZeros -o matrix_peaks_cut_AF488.gz --outFileSortedRegions regions_peaks_cut_AF488.bed

plotProfile -m matrix_peaks_cut_AF488.gz  -out TSS_profile_cut_AF488.pdf  --refPointLabel peak --regionsLabel ES EMS MS LS --outFileSortedRegions location_cut_AF488.bed --plotFileFormat pdf --perGroup
plotHeatmap -m matrix_peaks_cut_AF488.gz  -out TSS_heatmap_cut_AF488.pdf  --refPointLabel peak --regionsLabel ES EMS MS LS --outFileSortedRegions location_cut_AF488.bed --plotFileFormat pdf --perGroup \
--samplesLabel x39 x49 xw11 kong --heatmapHeight 35 --heatmapWidth 5 --colorMap OrRd


