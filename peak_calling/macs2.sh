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

threshold=0.01
type='p'
mkdir -p /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_${type}${threshold}
out="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_${type}${threshold}"

#ZH11-control-1
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam  -n ZH11_3_ES \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_3_ES_out

sort -k8,8nr ${out}/ZH11_3_ES_out/ZH11_3_ES_peaks.narrowPeak > ${out}/ZH11_3_ES_out/ZH11_3_ES_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam  -n ZH11_3_MS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_3_MS_out

sort -k8,8nr ${out}/ZH11_3_MS_out/ZH11_3_MS_peaks.narrowPeak > ${out}/ZH11_3_MS_out/ZH11_3_MS_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam -n ZH11_3_LS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_3_LS_out

sort -k8,8nr ${out}/ZH11_3_LS_out/ZH11_3_LS_peaks.narrowPeak > ${out}/ZH11_3_LS_out/ZH11_3_LS_peaks_sort.narrowPeak

#ZH11-control-2
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam  -n ZH11_4_ES \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_4_ES_out

sort -k8,8nr ${out}/ZH11_4_ES_out/ZH11_4_ES_peaks.narrowPeak > ${out}/ZH11_4_ES_out/ZH11_4_ES_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-2_S99_sort_rmdup_rmor_q30.bam  -n ZH11_4_MS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_4_MS_out

sort -k8,8nr ${out}/ZH11_4_MS_out/ZH11_4_MS_peaks.narrowPeak > ${out}/ZH11_4_MS_out/ZH11_4_MS_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-3_S99_sort_rmdup_rmor_q30.bam -n ZH11_4_LS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_4_LS_out

sort -k8,8nr ${out}/ZH11_4_LS_out/ZH11_4_LS_peaks.narrowPeak > ${out}/ZH11_4_LS_out/ZH11_4_LS_peaks_sort.narrowPeak

#ZH11-2
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-d-sup_S99_sort_rmdup_rmor_q30.bam  -n ZH11_2_G1 \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_2_G1_out

sort -k8,8nr ${out}/ZH11_2_G1_out/ZH11_2_G1_peaks.narrowPeak > ${out}/ZH11_2_G1_out/ZH11_2_G1_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/bams/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam  -n ZH11_2_ES \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_2_ES_out 

sort -k8,8nr ${out}/ZH11_2_ES_out/ZH11_2_ES_peaks.narrowPeak > ${out}/ZH11_2_ES_out/ZH11_2_ES_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam  -n ZH11_2_MS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_2_MS_out 

sort -k8,8nr ${out}/ZH11_2_MS_out/ZH11_2_MS_peaks.narrowPeak > ${out}/ZH11_2_MS_out/ZH11_2_MS_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam  -n ZH11_2_LS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_2_LS_out 

sort -k8,8nr ${out}/ZH11_2_LS_out/ZH11_2_LS_peaks.narrowPeak > ${out}/ZH11_2_LS_out/ZH11_2_LS_peaks_sort.narrowPeak

#ZH11-1
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam  -n ZH11_1_ES \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_1_ES_out

sort -k8,8nr ${out}/ZH11_1_ES_out/ZH11_1_ES_peaks.narrowPeak > ${out}/ZH11_1_ES_out/ZH11_1_ES_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam  -n ZH11_1_MS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_1_MS_out

sort -k8,8nr ${out}/ZH11_1_MS_out/ZH11_1_MS_peaks.narrowPeak > ${out}/ZH11_1_MS_out/ZH11_1_MS_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/ZH11_1/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam  -n ZH11_1_LS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/ZH11_1_LS_out

sort -k8,8nr ${out}/ZH11_1_LS_out/ZH11_1_LS_peaks.narrowPeak > ${out}/ZH11_1_LS_out/ZH11_1_LS_peaks_sort.narrowPeak

#NIP
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-6_sort_rmdup_rmor_q30.bam  -n NIP_ES \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/NIP_ES_out

sort -k8,8nr ${out}/NIP_ES_out/NIP_ES_peaks.narrowPeak > ${out}/NIP_ES_out/NIP_ES_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-4_sort_rmdup_rmor_q30.bam  -n NIP_MS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/NIP_MS_out

sort -k8,8nr ${out}/NIP_MS_out/NIP_MS_peaks.narrowPeak > ${out}/NIP_MS_out/NIP_MS_peaks_sort.narrowPeak

macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-5_sort_rmdup_rmor_q30.bam  -n NIP_LS \
        --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.0E8 -${type} ${threshold} --outdir ${out}/NIP_LS_out

sort -k8,8nr ${out}/NIP_LS_out/NIP_LS_peaks.narrowPeak > ${out}/NIP_LS_out/NIP_LS_peaks_sort.narrowPeak
