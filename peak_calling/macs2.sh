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

for file in $(find ./ -name '*q30.bam');do
        name=$(basename ${file})
        prefix=${name%%_*}
        macs2 callpeak -t  ${name} -n ${prefix} --shift 100 --extsize 200 --nomodel -B --SPMR -g 3.3E8  -q 0.01 --outdir ${out}/${prefix}_out 2> ${out}/${prefix}.macs2.log
        sort -k8,8nr ${out}/${prefix}_out/${prefix}_peaks.narrowPeak > ${out}/${prefix}_out/${prefix}_peaks_sort.narrowPeak
        sleep 1
done



	