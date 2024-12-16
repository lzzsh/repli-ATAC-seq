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

threshold=1e-5
type='p'
mkdir -p /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_${type}${threshold}
out="/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_${type}${threshold}"

#x39
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/x39_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/kong_sort_rmdup_rmor_q30.bam -n x39_cut \
        --nomode -f BAMPE --nomodel -B --SPMR -g 3.0E8 --keep-dup all -${type} ${threshold} --outdir ${out}/x39_cut_out 

sort -k8,8nr ${out}/x39_cut_out/x39_cut_peaks.narrowPeak > ${out}/x39_cut_out/x39_cut_peaks_sort.narrowPeak

#x49
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/x49_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/kong_sort_rmdup_rmor_q30.bam -n x49_cut \
        --nomode -f BAMPE --nomodel -B --SPMR -g 3.0E8 --keep-dup all -${type} ${threshold} --outdir ${out}/x49_cut_out 

sort -k8,8nr ${out}/x49_cut_out/x49_cut_peaks.narrowPeak > ${out}/x49_cut_out/x49_cut_peaks_sort.narrowPeak

#xw11
macs2 callpeak -t /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/xw11_sort_rmdup_rmor_q30.bam -c /storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/bwa_all_rawdata/kong_sort_rmdup_rmor_q30.bam -n xw11_cut \
        --nomode -f BAMPE --nomodel -B --SPMR -g 3.0E8 --keep-dup all -${type} ${threshold} --outdir ${out}/xw11_cut_out 

sort -k8,8nr ${out}/xw11_cut_out/xw11_cut_peaks.narrowPeak > ${out}/xw11_cut_out/xw11_cut_peaks_sort.narrowPeak

