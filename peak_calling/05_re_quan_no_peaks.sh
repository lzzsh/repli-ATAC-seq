#!/bin/bash

# Step 1: Extend peaks by 500bp
awk '{print $1, $2-500, $3+500}' OFS="\t" peaks_quan_tpm_filtered_pos_org.txt > peaks_extended.bed
awk '{if ($2 < 0) $2=0; print $0}' OFS="\t" peaks_extended.bed > peaks_extended_fixed.bed

# Step 2: Merge overlapping 500bp extended regions
bedtools sort -i peaks_extended_fixed.bed | bedtools merge -i - > peaks_extended_merged.bed

# Step 3: Compute read counts for 500bp extended regions
bedtools multicov -bams \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-d-sup_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-1_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-3_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-6_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-4_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-5_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-p_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-q_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-r_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/LZZ-2-g_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-h_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-i_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-j_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-k_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-l_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-2_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-3_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup2/bwa_all_rawdata/LZZ-10-S_S99_sort_rmdup_rmor_q30.bam \
-bed peaks_extended_merged.bed > peaks_open_region_counts.txt

# Step 4: Extend peaks by 5000bp
awk '{center=int(($2+$3)/2); print $1, center-5000, center+5000}' OFS="\t" peaks_quan_tpm_filtered_pos_org.txt > peaks_5kb_window.bed
awk '{if ($2 < 0) $2=0; print $0}' OFS="\t" peaks_5kb_window.bed > peaks_5kb_window_fixed.bed

# Step 5: Compute read counts for 5000bp extended regions
bedtools multicov -bams \
/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-d-sup_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-1-d_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-e_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-1-f_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-1_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-4_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-5_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-2_origin-2-6_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-3_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-6_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-4_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/01.RawData/origin-origin-5_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-p_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-q_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-r_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-2-g_result/bwa_all_rawdata/LZZ-2-g_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-h_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/bwa_all_rawdata/LZZ-2-i_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup/LZZ-adjp_result/merged_bam/LZZ-3-j_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-k_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_CR2/bwa_all_rawdata/LZZ-4-l_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-1_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-2_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-8-3_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-1_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-2_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata/LZZ-9-3_S99_sort_rmdup_rmor_q30.bam \
    /storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_sup2/bwa_all_rawdata/LZZ-10-S_S99_sort_rmdup_rmor_q30.bam \
-bed peaks_5kb_window_fixed.bed > peaks_5kb_window_counts.txt

sort -k1,1 -k2,2n peaks_5kb_window_counts.txt > peaks_5kb_window_counts_sorted.txt
sort -k1,1 -k2,2n peaks_open_region_counts.txt > peaks_open_region_counts_sorted.txt

# Step 6: Find overlaps between 5000bp regions and 500bp open regions
bedtools intersect -a peaks_5kb_window_fixed.bed -b peaks_extended_merged.bed -wo > peaks_5kb_vs_500bp_overlap.txt

# Step 7: Filter regions where overlap length is ≥ half of the open region length
awk '
{
    peak_length=$6-$5;
    if ($7 >= peak_length/2)  
        print $1, $2, $3, $4, $5, $6, $7;
}' OFS="\t" peaks_5kb_vs_500bp_overlap.txt > high_overlap_regions.txt

# Step 8: Extract counts of highly overlapping open regions
aawk '
BEGIN {
    # Read peaks_open_region_counts.txt and store counts
    while ((getline < "peaks_open_region_counts.txt") > 0) {
        key = $1"\t"$2"\t"$3;  # Create a key using chr, start, end
        open_counts[key] = "";  # Initialize an empty string for storing counts
        for (i=4; i<=NF; i++) {
            open_counts[key] = open_counts[key] "\t" $i;  # Concatenate all counts into a single string
        }
    }
}
{
    key = $4"\t"$5"\t"$6;  # Key for the high overlap region (500bp open region)
    
    if (key in open_counts) {
        print $1, $2, $3 open_counts[key];  # Directly output the 5000bp region with its associated counts
    } else {
        printf "%s\t%s\t%s", $1, $2, $3;  # Print the 5000bp window coordinates
        for (i=4; i<=NF; i++) {
            printf "\t0";  # If no matching counts are found, print 0
        }
        print "";  # End the line
    }
}' OFS="\t" high_overlap_regions.txt > high_overlap_counts.txt

# Step 9: Subtract overlapping open region counts from 5000bp counts
awk '
BEGIN {
    # Read high_overlap_counts.txt and store counts
    while ((getline < "high_overlap_counts.txt") > 0) {
        key = $1"\t"$2"\t"$3;  # Create key using chr, start, end
        overlap_counts[key] = "";  # Initialize an empty string
        for (i=4; i<=NF; i++) {
            overlap_counts[key] = overlap_counts[key] "\t" $i;  # Store entire count row as a string
        }
    }
}
{
    key = $1"\t"$2"\t"$3;  # Key for 5000bp window region

    printf "%s\t%s\t%s", $1, $2, $3;  # Print coordinates

    if (key in overlap_counts) {
        split(overlap_counts[key], overlap_values, "\t");  # Split stored values into an array
        for (i=4; i<=NF; i++) {
            adj_count = $i - overlap_values[i-3];  # Subtract corresponding value
            if (adj_count < 0) adj_count = 0;  # Ensure no negative values
            printf "\t%d", adj_count;
        }
    } else {
        for (i=4; i<=NF; i++) {
            printf "\t%d", $i;  # If no overlap, keep original count
        }
    }
    print "";  # End the line
}' OFS="\t" peaks_5kb_window_counts.txt > peaks_5kb_window_counts_filtered.txt

# Step 10: Compute final masked counts
paste peaks_5kb_window_counts_filtered.txt peaks_open_region_counts.txt | \
awk '{printf "%s\t%s\t%s", $1, $2, $3; for(i=4; i<=NF/2; i++) printf "\t%d", $i-$(i+NF/2); print ""}' > peaks_masked_counts.txt
