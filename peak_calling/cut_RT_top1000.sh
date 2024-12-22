sort -k8,8gr ../x39_cut_out/x39_cut_peaks.narrowPeak | head -n 2000 > ../x39_cut_out/x39_cut_peaks_top2000.narrowPeak
sort -k8,8gr ../x49_cut_out/x49_cut_peaks.narrowPeak | head -n 2000 > ../x49_cut_out/x49_cut_peaks_top2000.narrowPeak
sort -k8,8gr ../xw11_cut_out/xw11_cut_peaks.narrowPeak | head -n 2000 > ../xw11_cut_out/xw11_cut_peaks_top2000.narrowPeak

bedtools intersect -a ZH11_2_RT_control.gff3 -b ../x39_cut_out/x39_cut_peaks_top2000.narrowPeak -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $5, $6+$14, $6+$14+1, $4}' > x39_RT_top2000.bed
bedtools intersect -a ZH11_2_RT_control.gff3 -b ../x49_cut_out/x49_cut_peaks_top2000.narrowPeak -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $5, $6+$14, $6+$14+1, $4}' > x49_RT_top2000.bed
bedtools intersect -a ZH11_2_RT_control.gff3 -b ../xw11_cut_out/xw11_cut_peaks_top2000.narrowPeak -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $5, $6+$14, $6+$14+1, $4}' > xw11_RT_top2000.bed