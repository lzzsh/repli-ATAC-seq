# Set working directory
setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/peaks_quan/")

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
                     "ZH11-1-ES", "ZH11-1-MS", "ZH11-1-LS",
                     "NIP-1-ES", "NIP-1-MS", "NIP-1-LS")

# Function to normalize columns 4 to 19 using TPM
normalize_tpm <- function(data) {
  # Calculate region length
  data$RegionLength <- data$end - data$start + 1
  
  # Apply TPM normalization for each selected column
  data[, 4:19] <- lapply(data[, 4:19], function(col) {
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

# Apply normalization
normalized_count <- normalize_tpm(count)

# Save normalized data to a new file
output_file <- "./results/peaks_quan_tpm_normalized.txt"
write.table(normalized_count[,4:19], file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Print success message
cat("TPM normalization completed. Normalized data saved to:", output_file, "\n")