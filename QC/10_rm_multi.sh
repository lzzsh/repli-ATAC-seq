#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -J rm_multi
#SBATCH -o rm_multi.%j.out
#SBATCH -e rm_multi.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liaozizhuo@westlake.edu.cn

# Loop through all .bam files in the current directory
for file in *.bam; do
    echo "Processing $file ..."

    # Create a temporary file to store the filtered results
    temp_file="${file%.bam}_temp.bam"

    # Use samtools and grep to filter reads
    samtools view -h "$file" | \
    grep -v "chrUn" | \
    grep -v "chrSy" | \
    samtools view -bF 256 -o "$temp_file"

    # Replace the original file with the filtered results
    mv "$temp_file" "$file"

    echo "$file has been filtered and updated."
done

echo "All BAM files have been processed."
