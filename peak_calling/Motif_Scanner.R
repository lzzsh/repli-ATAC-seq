# Load Required R Packages**
library(TFBSTools)       # Handles transcription factor binding sites and motifs
library(JASPAR2024)      # Retrieves motifs from the JASPAR database
library(Biostrings)      # Handles DNA sequences
library(GenomicRanges)   # Handles genomic coordinates
library(Rsamtools)       # Reads indexed genome files
library(RSQLite)         # Connects to SQLite databases
library(tidyr)           # Data wrangling and reshaping
library(dplyr)           # Data manipulation

# Connect to JASPAR2024 and Retrieve the WOX11 Motif**
JASPAR2024 = JASPAR2024()
JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
query_result <- dbGetQuery(JASPARConnect, "SELECT * FROM MATRIX WHERE name LIKE '%WOX11%'")

# Select the first matching motif
if (nrow(query_result) > 0) {
  wox11_id <- query_result$ID[1]
  cat("Wox11 Motif ID:", wox11_id, "\n")
} else {
  stop("Error: No WOX11 motif found in JASPAR2024 database.")
}

# Retrieve PWM Data for the Selected Motif**
pfm_query <- dbGetQuery(JASPARConnect, paste0("SELECT * FROM MATRIX_DATA WHERE ID = '", wox11_id, "'"))

# Reformat the PWM Data
pfm_wide <- pfm_query %>%
  dplyr::select(row, col, val) %>%   # Select relevant columns
  pivot_wider(names_from = row, values_from = val) %>%  # Convert to wide format
  arrange(col)  # Ensure correct motif order

# Convert to matrix format and transpose for PWM creation
pfm_matrix <- as.matrix(pfm_wide[, -1])  # Remove `col` column
pfm_matrix <- t(round(pfm_matrix))   # Convert to integer counts

# Step 2: Define the pseudocount value
pseudocount <- 1  # You can adjust this value

# Step 3: Compute the modified Position Probability Matrix (PPM) using the given formula
col_sums <- colSums(pfm_matrix)  # Total count at each position
ppm_matrix <- sweep(pfm_matrix + pseudocount, 2, col_sums + 4 * pseudocount, "/")

# Step 4: Convert PPM to Position Weight Matrix (PWM)
background_freq <- c(0.25, 0.25, 0.25, 0.25)  # Assume uniform background
wox11_pwm <- log2(ppm_matrix / background_freq)  # Log-odds transformation

# Load Genome and Extract Sequences**
genome <- FaFile("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/rice_all_genomes_v7.fasta")

# Read peak regions from the narrowPeak file
peaks <- read.table("/storage/liuxiaodongLab/liaozizhuo/Projects//cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_cut_peaks_top1000.narrowPeak", header=FALSE)
gr <- GRanges(seqnames=peaks$V1, ranges=IRanges(start=peaks$V2, end=peaks$V3))

# Retrieve sequences corresponding to peaks
seqs <- getSeq(genome, gr)

# Perform PWM Motif Search**
# Ensure `seqs` is a DNAStringSet
seqs <- DNAStringSet(seqs)

# Scan sequences for motif occurrences
hits <- lapply(as.list(seqs), function(seq) matchPWM(wox11_pwm, seq, min.score = "80%"))

# Select peaks containing the motif
selected_indices <- which(lengths(hits) > 0)
selected_peaks <- peaks[selected_indices, ]

# Save Results**
write.table(selected_peaks, file="wox11_peaks_with_motif.bed", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("Total number of peaks containing the Wox11 motif:", nrow(selected_peaks), "\n")
