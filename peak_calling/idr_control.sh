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

cd /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/
out="/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/idr_result"

idr --samples macs2/NIP_ES_out/NIP_ES_peaks_sort.narrowPeak \
			  macs2/LZZ-1-d_out/LZZ-1-d_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/ES-idr --plot --log-output-file ${out}/ES.idr.log

idr --samples macs2/ZH11_1_MS_out/ZH11_1_MS_peaks_sort.narrowPeak \
			  macs2/LZZ-2-e_out/LZZ-2-e_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/MS-idr --plot --log-output-file ${out}/MS.idr.log

idr --samples macs2/ZH11_1_LS_out/ZH11_1_LS_peaks_sort.narrowPeak \
			  macs2/LZZ-1-f_out/LZZ-1-f_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/LS-idr --plot --log-output-file ${out}/LS.idr.log