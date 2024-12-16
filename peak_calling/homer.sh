# preparation
loadGenome.pl -name rice -fasta ./rice_all_genomes_v7.fasta -gtf ./all_DIY.gff3 -org "Oryza sativa" -version custom -force

samtools faidx rice_all_genomes_v7.fasta
cut -f1,2 rice_all_genomes_v7.fasta.fai > genome.chrom.sizes

sort -k8,8gr xw11_cut_peaks.narrowPeak | head -n 1000 | awk 'BEGIN{FS=OFS="\t"} {print $1, $2+$10, $2+$10+1, $4, $5, $6}' > xw11_summits_top1000.bed
bedtools slop -i xw11_summits_top1000.bed -g genome.chrom.sizes -b 150 > xw11_summits_top1000_temp.bed
mv xw11_summits_top1000_temp.bed xw11_summits_top1000.bed

# homer
findMotifsGenome.pl xw11_summits_top1000.bed rice ./output/