sort -k8,8gr ../x39_cut_out/x39_cut_peaks.narrowPeak | head -n 1000 > ../x39_cut_out/x39_cut_peaks_top1000.narrowPeak
sort -k8,8gr ../x49_cut_out/x49_cut_peaks.narrowPeak | head -n 1000 > ../x49_cut_out/x49_cut_peaks_top1000.narrowPeak
sort -k8,8gr ../xw11_cut_out/xw11_cut_peaks.narrowPeak | head -n 1000 > ../xw11_cut_out/xw11_cut_peaks_top1000.narrowPeak

# cuttag site
bedtools intersect -a ZH11_RT_all.gff3 -b ../x39_cut_out/x39_cut_peaks_top1000.narrowPeak -wb | awk 'BEGIN{FS=OFS="\t"} {print $5, $6+$14, $6+$14+1, $4}' > ../x39_cut_out/x39_RT_top1000.bed
bedtools intersect -a ZH11_RT_all.gff3 -b ../x49_cut_out/x49_cut_peaks_top1000.narrowPeak  -wb | awk 'BEGIN{FS=OFS="\t"} {print $5, $6+$14, $6+$14+1, $4}' > ../x49_cut_out/x49_RT_top1000.bed
bedtools intersect -a ZH11_RT_all.gff3 -b ../xw11_cut_out/xw11_cut_peaks_top1000.narrowPeak -wb | awk 'BEGIN{FS=OFS="\t"} {print $5, $6+$14, $6+$14+1, $4}' > ../xw11_cut_out/xw11_RT_top1000.bed

# 
bedtools intersect -a ZH11_RT.gff3 -b ../x39_cut_out/x39_cut_peaks_top1000.narrowPeak -wa -wb > ../x39_cut_out/x39_RT_top1000.bed
bedtools intersect -a ZH11_RT.gff3 -b ../x49_cut_out/x49_cut_peaks_top1000.narrowPeak -wa -wb > ../x49_cut_out/x49_RT_top1000.bed
bedtools intersect -a ZH11_RT.gff3 -b ../xw11_cut_out/xw11_cut_peaks_top1000.narrowPeak -wa -wb > ../xw11_cut_out/xw11_RT_top1000.bed

bedtools intersect -a ZH11_RT.gff3 -b ../x39_cut_out/x39_cut_peaks.narrowPeak -wa -wb > ../x39_cut_out/x39_RT.bed
bedtools intersect -a ZH11_RT.gff3 -b ../x49_cut_out/x49_cut_peaks.narrowPeak -wa -wb > ../x49_cut_out/x49_RT.bed
bedtools intersect -a ZH11_RT.gff3 -b ../xw11_cut_out/xw11_cut_peaks.narrowPeak -wa -wb > ../xw11_cut_out/xw11_RT.bed