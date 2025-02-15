# Set working directory
setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load data
count <- read.table("./peaks_quan.txt")

# Add column names
colnames(count) <- c("chr", "start", "end",
                     "ZH11-3-ES", "ZH11-3-MS", "ZH11-3-LS",
                     "ZH11-4-ES", "ZH11-4-MS", "ZH11-4-LS",
                     "ZH11-2-G1", "ZH11-2-ES", "ZH11-2-MS", "ZH11-2-LS",
                     "ZH11-1-G1", "ZH11-1-ES", "ZH11-1-MS", "ZH11-1-LS",
                     "NIP-1-G1", "NIP-1-ES", "NIP-1-MS", "NIP-1-LS")

# Directory to save the results
output_dir <- "./results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Remove the G1 column
count <- count %>% select(-c("ZH11-2-G1", "ZH11-1-G1", "NIP-1-G1"))

# Initialize data frames for plotting
all_data_no_outliers <- data.frame()
all_data_with_outliers_removed <- data.frame()

# Function to normalize columns 4 to 19 using TPM
normalize_tpm <- function(data) {
  # Calculate region length
  data$RegionLength <- data$end - data$start + 1
  
  # Apply TPM normalization for each selected column
  data[, 4:18] <- lapply(data[, 4:18], function(col) {
    # Calculate RPK (Reads Per Kilobase)
    rpk <- col / (data$RegionLength / 1000)
    # Calculate scaling factor (Total RPK / 10^6)
    scaling_factor <- sum(rpk) / 1e6
    # Calculate TPM
    tpm <- rpk / scaling_factor
    return(tpm)
  })
  
  return(data)
}

normalized_count <- normalized_count[rowMeans(normalized_count[4:18], na.rm = TRUE) < 100,]
org <- read.table("./organelle_peaks.bed")
org <- org %>% unique() %>% setNames(c("chr","start","end"))

# Remove rows from normalized_count that match the first three columns in org
normalized_count <- anti_join(normalized_count, org[, 1:3], by = colnames(org)[1:3])
