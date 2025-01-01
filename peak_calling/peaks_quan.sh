# sort and merge peaks
cat ../ZH11_3_ES_out/ZH11_3_ES_peaks.narrowPeak ../ZH11_3_MS_out/ZH11_3_MS_peaks_sort.narrowPeak ../ZH11_3_LS_out/ZH11_3_LS_peaks_sort.narrowPeak > ZH11_3_peaks_p0.01
cat ../ZH11_2_ES_out/ZH11_2_ES_peaks.narrowPeak ../ZH11_2_MS_out/ZH11_2_MS_peaks_sort.narrowPeak ../ZH11_2_LS_out/ZH11_2_LS_peaks_sort.narrowPeak > ZH11_2_peaks_p0.01
cat ZH11_3_peaks_p0.01 ZH11_2_peaks_p0.01 > peaks_p0.01
sort -k1,1 -k2,2n peaks_p0.01 > peaks_p0.01_sort
bedtools merge -i peaks_p0.01_sort > merged_peaks_p0.01_sort

# select reads only from open chromatin region
files=(
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-d-sup_S99_sort_rmdup_rmor_q30.bam"
	)

for file in ${files[@]}; do
	name=$(basename ${file})
	name=${name%%.*}
	bedtools intersect -abam ${file} -b merged_peaks_p0.01_sort -wa > ${name}_filtered.bam
    samtools index ${name}_filtered.bam
done

#mkdir bws
for file in $(find ./ -name '*_filtered.bam');do
	name=$(basename ${file})
	prefix=${name%_*}
	bamCoverage --bam ${file} -o ./bws/${prefix}.bw --normalizeUsing BPM -p 5 --binSize 50
	sleep 1
done