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

idr --samples ZH11_1/macs2/origin-2_origin-2-1_sort_rmdup_rmor_out/origin-2_origin-2-1_sort_rmdup_rmor_peaks_sort.narrowPeak \
			  macs2/LZZ-1-d-sup_out/LZZ-1-d-sup_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/G1-idr --plot --log-output-file ${out}/G1.idr.log

idr --samples ZH11_1/macs2/origin-2_origin-2-2_sort_rmdup_rmor_out/origin-2_origin-2-2_sort_rmdup_rmor_peaks_sort.narrowPeak \
			  macs2/LZZ-1-e-sup_out/LZZ-1-e-sup_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/MS_withouEdU-idr --plot --log-output-file ${out}/MS_withouEdU.idr.log

idr --samples ZH11_1/macs2/origin-2_origin-2-3_sort_rmdup_rmor_out/origin-2_origin-2-3_sort_rmdup_rmor_peaks_sort.narrowPeak \
			  macs2/LZZ-1-f-sup_out/LZZ-1-f-sup_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/G2-idr --plot --log-output-file ${out}/G2.idr.log

idr --samples ZH11_1/macs2/origin-2_origin-2-4_sort_rmdup_rmor_out/origin-2_origin-2-4_sort_rmdup_rmor_peaks_sort.narrowPeak \
			  macs2/LZZ-1-d_out/LZZ-1-d_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/ES-idr --plot --log-output-file ${out}/ES.idr.log

idr --samples ZH11_1/macs2/origin-2_origin-2-5_sort_rmdup_rmor_out/origin-2_origin-2-5_sort_rmdup_rmor_peaks_sort.narrowPeak \
			  macs2/LZZ-2-e_out/LZZ-2-e_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/MS-idr --plot --log-output-file ${out}/MS.idr.log

idr --samples ZH11_1/macs2/origin-2_origin-2-6_sort_rmdup_rmor_out/origin-2_origin-2-6_sort_rmdup_rmor_peaks_sort.narrowPeak \
			  macs2/LZZ-1-f_out/LZZ-1-f_peaks_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${out}/LS-idr --plot --log-output-file ${out}/LS.idr.log