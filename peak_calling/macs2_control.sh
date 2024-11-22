#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J macs2
#SBATCH -o macs2.%j.out
#SBATCH -e macs2.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata
mkdir ../macs2
out="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2"

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam -n LZZ-1-d \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -q 0.01 --outdir ${out}/LZZ-1-d_out 

sort -k8,8nr ${out}/LZZ-1-d_out/LZZ-1-d_peaks.narrowPeak > ${out}/LZZ-1-d_out/LZZ-1-d_peaks_sort.narrowPeak

macs2 callpeak -t LZZ-2-e_S99_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam -n LZZ-2-e \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -q 0.01 --outdir ${out}/LZZ-2-e_out 

sort -k8,8nr ${out}/LZZ-2-e_out/LZZ-2-e_peaks.narrowPeak > ${out}/LZZ-2-e_out/LZZ-2-e_peaks_sort.narrowPeak

macs2 callpeak -t LZZ-1-f_S99_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam -n LZZ-1-f \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -q 0.01 --outdir ${out}/LZZ-1-f_out 

sort -k8,8nr ${out}/LZZ-1-f_out/LZZ-1-f_peaks.narrowPeak > ${out}/LZZ-1-f_out/LZZ-1-f_peaks_sort.narrowPeak

