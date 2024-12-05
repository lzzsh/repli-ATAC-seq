#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J macs2
#SBATCH -o macs2.%j.out
#SBATCH -e macs2.%j.err
#SBATCH -p lxd_20241204,intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

threshold=0.01
type='q'
mkdir -p /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/macs2/macs2_${type}${threshold}
out="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/macs2/macs2_${type}${threshold}"

# wox11
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/LZZ-2-g_S99_sort_rmdup_rmor_q30.bam  -n xw11_CR_ES \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/xw11_CR_ES_out

sort -k8,8nr ${out}/xw11_CR_ES_out/xw11_CR_ES_peaks.narrowPeak > ${out}/xw11_CR_ES_out/xw11_CR_ES_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-h_S99_sort_rmdup_rmor_q30.bam  -n xw11_CR_MS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/xw11_CR_MS_out

sort -k8,8nr ${out}/xw11_CR_MS_out/xw11_CR_MS_peaks.narrowPeak > ${out}/xw11_CR_MS_out/xw11_CR_MS_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-i_S99_sort_rmdup_rmor_q30.bam  -n xw11_CR_LS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/xw11_CR_LS_out

sort -k8,8nr ${out}/xw11_CR_LS_out/xw11_CR_LS_peaks.narrowPeak > ${out}/xw11_CR_LS_out/xw11_CR_LS_peaks_sort.narrowPeak

# x39
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/LZZ-3-j_S99_sort_rmdup_rmor_q30.bam  -n x39_CR_ES \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/x39_CR_ES_out

sort -k8,8nr ${out}/x39_CR_ES_out/x39_CR_ES_peaks.narrowPeak > ${out}/x39_CR_ES_out/x39_CR_ES_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-k_S99_sort_rmdup_rmor_q30.bam  -n x39_CR_MS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/x39_CR_MS_out

sort -k8,8nr ${out}/x39_CR_MS_out/x39_CR_MS_peaks.narrowPeak > ${out}/x39_CR_MS_out/x39_CR_MS_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-l_S99_sort_rmdup_rmor_q30.bam  -n x39_CR_LS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/x39_CR_LS_out

sort -k8,8nr ${out}/x39_CR_LS_out/x39_CR_LS_peaks.narrowPeak > ${out}/x39_CR_LS_out/x39_CR_LS_peaks_sort.narrowPeak
