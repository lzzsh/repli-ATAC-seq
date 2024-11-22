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

macs2 callpeak -t ../ZH11-1/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam -n ZH11_1_ES \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -q 0.01 --outdir ${out}/ZH11_1_ES_out 

sort -k8,8nr ${out}/ZH11_1_ES_out/ZH11_1_ES_peaks.narrowPeak > ${out}/ZH11_1_ES_out/ZH11_1_ES_peaks_sort.narrowPeak

macs2 callpeak -t ../ZH11-1/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam -n ZH11_1_MS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -q 0.01 --outdir ${out}/ZH11_1_MS_out 

sort -k8,8nr ${out}/ZH11_1_MS_out/ZH11_1_MS_peaks.narrowPeak > ${out}/ZH11_1_MS_out/ZH11_1_MS_peaks_sort.narrowPeak

macs2 callpeak -t ../ZH11-1/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam -n ZH11_1_LS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -q 0.01 --outdir ${out}/ZH11_1_LS_out 

sort -k8,8nr ${out}/ZH11_1_LS_out/ZH11_1_LS_peaks.narrowPeak > ${out}/ZH11_1_LS_out/ZH11_1_LS_peaks_sort.narrowPeak

macs2 callpeak -t ../01.RawData/origin-origin-6_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam -n NIP_ES \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -q 0.01 --outdir ${out}/NIP_ES_out 

sort -k8,8nr ${out}/NIP_ES_out/NIP_ES_peaks.narrowPeak > ${out}/NIP_ES_out/NIP_ES_peaks_sort.narrowPeak