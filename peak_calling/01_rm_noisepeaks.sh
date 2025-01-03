#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J rm_noise
#SBATCH -o rm_noise.%j.out
#SBATCH -e rm_noise.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

module load bedtools/2.30.0
module load samtools/1.13

# calculate reads number in each open chromatin region
ZH11_1_location="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData"
ZH11_2_location="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata"
ZH11_3_location="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata"
ZH11_4_location="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata"
NIP_location="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData"
threshold=0.01
type='p'
macs="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_${type}${threshold}"

out="${macs}/rm_noisepeak"
mkdir ${out}

# ZH11-1
bedtools multicov -bams  ${ZH11_1_location}/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_1_ES_out/ZH11_1_ES_peaks_sort.narrowPeak > ZH11_1_ES.txt
bedtools multicov -bams  ${ZH11_1_location}/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_1_MS_out/ZH11_1_MS_peaks_sort.narrowPeak > ZH11_1_MS.txt
bedtools multicov -bams  ${ZH11_1_location}/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_1_LS_out/ZH11_1_LS_peaks_sort.narrowPeak > ZH11_1_LS.txt

# ZH11-2
bedtools multicov -bams  /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_2_ES_out/ZH11_1_ES_peaks_sort.narrowPeak > ZH11_2_ES.txt
bedtools multicov -bams  ${ZH11_2_location}/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_2_MS_out/ZH11_2_MS_peaks_sort.narrowPeak > ZH11_2_MS.txt
bedtools multicov -bams  ${ZH11_2_location}/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_2_LS_out/ZH11_2_LS_peaks_sort.narrowPeak > ZH11_2_LS.txt

# NIP
bedtools multicov -bams  ${NIP_location}/origin-origin-6_sort_rmdup_rmor_q30.bam -bed ${macs}/NIP_ES_out/NIP_ES_peaks_sort.narrowPeak > NIP_ES.txt
bedtools multicov -bams  ${NIP_location}/origin-origin-4_sort_rmdup_rmor_q30.bam -bed ${macs}/NIP_MS_out/NIP_MS_peaks_sort.narrowPeak > NIP_MS.txt
bedtools multicov -bams  ${NIP_location}/origin-origin-5_sort_rmdup_rmor_q30.bam -bed ${macs}/NIP_LS_out/NIP_LS_peaks_sort.narrowPeak > NIP_LS.txt

# ZH11-3(control-1)
bedtools multicov -bams  ${ZH11_3_location}/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_3_ES_out/ZH11_3_ES_peaks_sort.narrowPeak > ZH11_3_ES.txt
bedtools multicov -bams  ${ZH11_3_location}/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_3_MS_out/ZH11_3_MS_peaks_sort.narrowPeak > ZH11_3_MS.txt
bedtools multicov -bams  ${ZH11_3_location}/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_3_LS_out/ZH11_3_LS_peaks_sort.narrowPeak > ZH11_3_LS.txt

# ZH11-4(control-2)
bedtools multicov -bams  ${ZH11_4_location}/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_4_ES_out/ZH11_4_ES_peaks_sort.narrowPeak > ZH11_4_ES.txt
bedtools multicov -bams  ${ZH11_4_location}/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_4_MS_out/ZH11_4_MS_peaks_sort.narrowPeak > ZH11_4_MS.txt
bedtools multicov -bams  ${ZH11_4_location}/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam -bed ${macs}/ZH11_4_LS_out/ZH11_4_LS_peaks_sort.narrowPeak > ZH11_4_LS.txt