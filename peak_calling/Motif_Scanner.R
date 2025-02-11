# Load Required Packages
library(TFBSTools)      
library(JASPAR2024)     
library(Biostrings)      
library(GenomicRanges)   
library(Rsamtools)       
library(RSQLite)        
library(tidyr)           
library(dplyr)           

# Connect to JASPAR2024 and Retrieve WOX11 Motifs
JASPAR2024 <- JASPAR2024()
JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

# Query for motifs containing 'WOX11'
query_result <- dbGetQuery(JASPARConnect, "SELECT * FROM MATRIX WHERE name LIKE '%WOX11%'")

# Ensure at least two motifs are found
if (nrow(query_result) < 2) {
  stop("Error: Less than two WOX11 motifs found in JASPAR2024 database.")
}

# Extract the first two motif IDs
wox11_id1 <- query_result$ID[1]
wox11_id2 <- query_result$ID[2]

cat("WOX11 Motif IDs:", wox11_id1, "and", wox11_id2, "\n")

# Function to retrieve, process, and compute PWM for a given motif ID
get_pwm <- function(motif_id) {
  pfm_query <- dbGetQuery(JASPARConnect, paste0("SELECT * FROM MATRIX_DATA WHERE ID = '", motif_id, "'"))
  
  pfm_wide <- pfm_query %>%
    dplyr::select(row, col, val) %>%
    pivot_wider(names_from = row, values_from = val) %>%
    arrange(col)
  
  pfm_matrix <- as.matrix(pfm_wide[, -1])
  pfm_matrix <- t(round(pfm_matrix))
  
  # Apply Pseudocounts and Compute PWM
  pseudocount <- 1  
  col_sums <- colSums(pfm_matrix)  
  ppm_matrix <- sweep(pfm_matrix + pseudocount, 2, col_sums + 4 * pseudocount, "/")
  
  background_freq <- c(0.25, 0.25, 0.25, 0.25)  
  pwm_matrix <- log2(ppm_matrix / background_freq)  
  
  return(pwm_matrix)
}

# Compute PWM for both motifs
wox11_pwm1 <- get_pwm(wox11_id1)
wox11_pwm2 <- get_pwm(wox11_id2)

# Load Genome and Extract Sequences
genome <- FaFile("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta")

# Read peak regions from the narrowPeak file
peaks <- read.table("/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_cut_peaks_top1000.narrowPeak", header=FALSE)

# Assign column names for clarity (assuming narrowPeak format)
colnames(peaks) <- c("chr", "start", "end", "V4", "score", "strand", "signalValue", "pValue", "qValue", "peak")

# Convert to GRanges object
gr <- GRanges(seqnames=peaks$chr, ranges=IRanges(start=peaks$start, end=peaks$end))

# Retrieve sequences corresponding to peaks
seqs <- getSeq(genome, gr)
seqs <- DNAStringSet(seqs)  # Ensure sequences are in DNAStringSet format

# Function to scan sequences for a given PWM
scan_pwm_hits <- function(pwm, seqs, threshold="70%") {
  hits <- lapply(as.list(seqs), function(seq) matchPWM(pwm, seq, min.score = threshold))
  selected_indices <- which(lengths(hits) > 0)
  return(selected_indices)
}

# Scan sequences for both motifs
selected_indices1 <- scan_pwm_hits(wox11_pwm1, seqs, threshold="80%")
selected_indices2 <- scan_pwm_hits(wox11_pwm2, seqs, threshold="80%")

# Extract peaks containing each motif
selected_peaks1 <- peaks[selected_indices1, ]
selected_peaks2 <- peaks[selected_indices2, ]

# Save motif-specific results
write.table(selected_peaks1, file="/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/wox11_peaks_with_motif1.bed", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(selected_peaks2, file="/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/wox11_peaks_with_motif2.bed", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# Merge results and remove duplicates based on V4 (peak ID)
combined_unique_peaks <- rbind(selected_peaks1, selected_peaks2) %>%
  distinct(V4, .keep_all = TRUE)  # Ensure V4 is unique while keeping all columns

# Save the merged unique peaks
write.table(combined_unique_peaks, file="/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/combined_unique_peaks.bed", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# Print Summary
cat("Total peaks containing Motif 1:", nrow(selected_peaks1), "\n")
cat("Total peaks containing Motif 2:", nrow(selected_peaks2), "\n")
cat("Total unique peaks (after merging):", nrow(combined_unique_peaks), "\n")