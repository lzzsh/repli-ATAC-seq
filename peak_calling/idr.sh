#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J idr
#SBATCH -o idr.%j.out
#SBATCH -e idr.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

threshold=0.1
type='p'
mkdir /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/idr_result
input="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_${type}${threshold}"
out="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/idr_result"

# ZH11-1 and ZH11-2
# idr --samples ${input}/NIP_ES_out/NIP_ES_peaks_sort.narrowPeak \
# 			  ${input}/LZZ-1-d_out/LZZ-1-d_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/ES-idr --plot --log-output-file ${out}/ES.idr.log

# idr --samples ${input}/ZH11_1_MS_out/ZH11_1_MS_peaks_sort.narrowPeak \
# 			  ${input}/LZZ-2-e_out/LZZ-2-e_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/MS-idr --plot --log-output-file ${out}/MS.idr.log

# idr --samples ${input}/ZH11_1_LS_out/ZH11_1_LS_peaks_sort.narrowPeak \
# 			  ${input}/LZZ-1-f_out/LZZ-1-f_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/LS-idr --plot --log-output-file ${out}/LS.idr.log

# ZH11-1 and NIP
idr --samples ${input}/ZH11_1_ES_out/ZH11_1_ES_peaks_sort.narrowPeak \
                          ${input}/NIP_ES_out/NIP_ES_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/ES-idr --plot --log-output-file ${out}/ES.idr.log

idr --samples ${input}/ZH11_1_MS_out/ZH11_1_MS_peaks_sort.narrowPeak \
                          ${input}/NIP_MS_out/NIP_MS_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/MS-idr --plot --log-output-file ${out}/MS.idr.log

idr --samples ${input}/ZH11_1_LS_out/ZH11_1_LS_peaks_sort.narrowPeak \
                          ${input}/NIP_LS_out/NIP_LS_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/LS-idr --plot --log-output-file ${out}/LS.idr.log
