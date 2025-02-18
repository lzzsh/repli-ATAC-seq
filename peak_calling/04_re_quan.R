# Set working directory
setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load data
count <- read.table("./peaks_re_quan.txt")

# Add column names
colnames(count) <- c("chr", "start", "end",
                     "ZH11-3-ES", "ZH11-3-MS", "ZH11-3-LS",
                     "ZH11-4-ES", "ZH11-4-MS", "ZH11-4-LS",
                     "ZH11-2-G1", "ZH11-2-ES", "ZH11-2-MS", "ZH11-2-LS",
                     "ZH11-1-G1", "ZH11-1-ES", "ZH11-1-MS", "ZH11-1-LS",
                     "NIP-1-G1", "NIP-1-ES", "NIP-1-MS", "NIP-1-LS")

# Remove the G1 column
count <- count %>% select(-c("ZH11-2-G1", "ZH11-1-G1", "NIP-1-G1"))

# Directory to save the results
output_dir <- "."
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to normalize columns 4 to 19 using TPM
normalize_tpm <- function(data) {
  # Calculate region length
  data$RegionLength <- data$end - data$start + 1
  
  # Apply TPM normalization for each selected column
  data[, 4:18] <- lapply(data[, 4:18], function(col) {
    # Calculate RPK (Reads Per Kilob base)
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

# Prepare data for plotting
all_data <- do.call(rbind, lapply(4:18, function(i) {
  data.frame(value = log(normalized_count[, i] + 1), sample = colnames(normalized_count)[i])
}))

# Save normalized data to a new file
output_file <- "./peaks_re_quan_tpm_normalized.txt"
write.table(normalized_count[,1:19], file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Create a boxplot for all samples
p <- ggplot(all_data, aes(x = sample, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12)) +
  ggtitle("Boxplots for All Samples (Normalized Data)") +
  ylab("Log(TPM + 1)") +
  xlab("Sample")
ggsave(file.path(output_dir, "boxplots_normalized_combined_requan.pdf"), p, width = 12, height = 6)

# Print success message
cat("TPM normalization and visualization completed.\n")
cat("Normalized data saved to:", output_file, "\n")
