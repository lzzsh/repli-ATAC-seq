# sort and merge peaks
cat ../ZH11_1_ES_out/ZH11_1_ES_peaks.narrowPeak ../ZH11_1_MS_out/ZH11_1_MS_peaks_sort.narrowPeak ../ZH11_1_LS_out/ZH11_1_LS_peaks_sort.narrowPeak > ZH11_1_peaks_p0.05
cat ../ZH11_2_ES_out/ZH11_2_ES_peaks.narrowPeak ../ZH11_2_MS_out/ZH11_2_MS_peaks_sort.narrowPeak ../ZH11_2_LS_out/ZH11_2_LS_peaks_sort.narrowPeak > ZH11_2_peaks_p0.05
cat ZH11_1_peaks_p0.05 ZH11_2_peaks_p0.05 > peaks_p0.05
sort -k1,1 -k2,2n peaks_p0.05 > peaks_p0.05_sort
bedtools merge -i peaks_p0.05_sort > merged_peaks_p0.05_sort

# select reads only from open chromatin region
files=(
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam"
	)

for file in ${files[@]}; do
	name=$(basename ${file})
	name=${name%%.*}
	bedtools intersect -abam ${file} -b merged_peaks_p0.05_sort -wa > ${name}_filtered.bam
    samtools index ${name}_filtered.bam
done

#mkdir bws
for file in $(find ./ -name '*_filtered.bam');do
	name=$(basename ${file})
	prefix=${name%_*}
	bamCoverage --bam ${file} -o ./bws/${prefix}.bw --normalizeUsing BPM -p 5 --binSize 50
	sleep 1
done