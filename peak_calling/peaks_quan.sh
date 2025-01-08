# sort and merge peaks
awk '{print $1"\t"$2"\t"$3}' *.txt > all_peaks.txt
sort -k1,1 -k2,2n all_peaks.txt | bedtools merge -i - > merged_peaks.bed

# select reads only from open chromatin region
files=(
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-2_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-3_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-d-sup_S99_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-4_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-5_sort_rmdup_rmor_q30.bam"
	"/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-6_sort_rmdup_rmor_q30.bam"
	)

for file in ${files[@]}; do
	name=$(basename ${file})
	name=${name%%.*}
	bedtools intersect -abam ${file} -b merged_peaks.bed -wa > ${name}_filtered.bam
    samtools index ${name}_filtered.bam
done

#mkdir bws
for file in $(find ./ -name '*_filtered.bam');do
	name=$(basename ${file})
	prefix=${name%_*}
	bamCoverage --bam ${file} -o ./bws/${prefix}.bw --normalizeUsing BPM -p 5 --binSize 50
	sleep 1
done