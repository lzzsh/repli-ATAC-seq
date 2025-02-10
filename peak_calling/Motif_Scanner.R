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
pwm_query <- dbGetQuery(JASPARConnect, paste0("SELECT * FROM MATRIX_DATA WHERE ID = '", wox11_id, "'"))

# Reformat the PWM Data
pwm_wide <- pwm_query %>%
  dplyr::select(row, col, val) %>%   # Select relevant columns
  pivot_wider(names_from = row, values_from = val) %>%  # Convert to wide format
  arrange(col)  # Ensure correct motif order

# Convert to matrix format and transpose for PWM creation
pwm_matrix <- as.matrix(pwm_wide[, -1])  # Remove `col` column
pfm_matrix <- round(pwm_matrix * 1000)   # Convert to integer counts
storage.mode(pfm_matrix) <- "integer"    # Ensure matrix is integer
rownames(pfm_matrix) <- c("A", "C", "G", "T")  # Assign row names

# Convert PFM to PWM
wox11_pwm <- PWM(t(pfm_matrix), type = "log2probratio",
                 prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25))

# Load Genome and Extract Sequences**
genome <- FaFile("/Users/lzz/Desktop/repli-ATAC-seq/reference_genome/rice_all_genomes_v7.fasta")

# Read peak regions from the narrowPeak file
peaks <- read.table("/Users/lzz/Desktop/repli-ATAC-seq/cuttag/macs2/xw11_cut_peaks_top1000.narrowPeak", header=FALSE)
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