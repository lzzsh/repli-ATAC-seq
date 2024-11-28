# select the peaks with IDR < 0.05
awk '{if($5 >= 540) print $0}' ES-idr > ES-idr_0.05
awk '{if($5 >= 540) print $0}' MS-idr > MS-idr_0.05
awk '{if($5 >= 540) print $0}' LS-idr > LS-idr_0.05

# sort and merge peaks
cat ES-idr_0.05 MS-idr_0.05 LS-idr_0.05 > EdU-idr_0.05
sort -k1,1 -k2,2n EdU-idr_0.05 > EdU-idr_0.05_sort
bedtools merge -i EdU-idr_0.05_sort > EdU-idr_0.05_sort_merge

# calculate reads number in each open chromatin region
location="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData"
bedtools multicov -bams \
${location}/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam \
${location}/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam \
${location}/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam \
-bed EdU-idr_0.05_sort_merge > EdU-idr_reads_ZH11.txt

# select the peaks with IDR < 0.1
awk '{if($5 >= 415) print $0}' ES-idr > ES-idr_0.1
awk '{if($5 >= 415) print $0}' MS-idr > MS-idr_0.1
awk '{if($5 >= 415) print $0}' LS-idr > LS-idr_0.1

# sort and merge peaks
cat ES-idr_0.1 MS-idr_0.1 LS-idr_0.1 > EdU-idr_0.1
sort -k1,1 -k2,2n EdU-idr_0.1 > EdU-idr_0.1_sort
bedtools merge -i EdU-idr_0.1_sort > EdU-idr_0.1_sort_merge

# calculate reads number in each open chromatin region
location="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData"
bedtools multicov -bams \
${location}/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam \
${location}/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam \
${location}/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam \
-bed EdU-idr_0.1_sort_merge > EdU-idr_reads_ZH11.txt