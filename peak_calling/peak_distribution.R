# Load required R packages
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

# Step 1: Read the GFF3 file and extract gene annotation information
gff_file <- "/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3"
gff_data <- import(gff_file, format="gff3")

# Step 2: Read WOX11 CUT&Tag peak data
peaks_file <- "/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_cut_peaks_top1000.narrowPeak"
peaks <- read.table(peaks_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(peaks) <- c("chr", "start", "end", "V4", "score", "strand", "signalValue", "pValue", "qValue", "peak")

# Convert peaks to GenomicRanges
gr_peaks <- GRanges(seqnames = peaks$chr, 
                    ranges = IRanges(start = peaks$start, end = peaks$end))

# Step 3: Extract different genomic regions from GFF3 annotations
# Extract genes, mRNAs, 5' UTRs, and 3' UTRs
genes <- gff_data[gff_data$type == "gene"]
mrnas <- gff_data[gff_data$type == "mRNA"]
utr5 <- gff_data[gff_data$type == "five_prime_UTR"]
utr3 <- gff_data[gff_data$type == "three_prime_UTR"]

# Step 4: Define promoter regions (2000 bp upstream of gene start sites)
promoters <- genes
start(promoters) <- start(promoters) - 2000  # Define promoter upstream region
end(promoters) <- start(promoters) + 2000   # Define promoter region spanning 4000 bp

# Step 5: Compute overlaps between peaks and different genomic regions
overlap_promoter <- countOverlaps(gr_peaks, promoters) > 0
overlap_utr5 <- countOverlaps(gr_peaks, utr5) > 0
overlap_utr3 <- countOverlaps(gr_peaks, utr3) > 0
overlap_gene <- countOverlaps(gr_peaks, genes) > 0

# Identify intergenic peaks (not overlapping any known gene-related regions)
overlap_intergenic <- !(overlap_promoter | overlap_utr5 | overlap_utr3 | overlap_gene)

# Step 6: Count the number of peaks in each genomic category
peak_counts <- data.frame(
  Category = c("Promoter", "5' UTR", "3' UTR", "Gene", "Intergenic"),
  Count = c(sum(overlap_promoter), sum(overlap_utr5), sum(overlap_utr3), sum(overlap_gene), sum(overlap_intergenic))
)

# Calculate the percentage of peaks in each region
peak_counts <- peak_counts %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Step 7: Display the peak distribution results
print(peak_counts)

# Step 8: Generate a bar plot for peak distribution using Morandi color scheme
ggplot(peak_counts, aes(x=Category, y=Percentage, fill=Category)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#D4E4C3", "#EAC9C1", "#F0DBCD", "#D6C4C4", "#E8E6C4")) +  # Apply Morandi colors
  theme_minimal() +
  labs(title="WOX11 CUT&Tag Peaks Distribution", y="Percentage (%)", x="Genomic Region") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme_classic()