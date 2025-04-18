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

# Function to remove local outliers (within each sample) using IQR
remove_local_outliers <- function(data) {
  signal_cols <- 4:18
  initial_rows <- nrow(data)
  
  # Data frame to store outlier information
  outlier_info <- data.frame()
  
  # Apply IQR-based filtering to each signal column
  for (col_idx in signal_cols) {
    col_name <- colnames(data)[col_idx]
    Q1 <- quantile(data[[col_idx]], 0.25, na.rm = TRUE)
    Q3 <- quantile(data[[col_idx]], 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 2.0 * IQR
    upper_bound <- Q3 + 2.0 * IQR
    
    # Identify outliers
    outliers <- (data[[col_idx]] < lower_bound) | (data[[col_idx]] > upper_bound)
    
    # Store outlier information
    if (any(outliers, na.rm = TRUE)) {
      outlier_data <- data[outliers, c("chr", "start", "end")]
      outlier_data$Sample <- col_name
      outlier_data$OutlierValue <- data[outliers, col_idx]
      outlier_info <- rbind(outlier_info, outlier_data)
    }
    
    # Set outliers to NA
    data[[col_idx]][outliers] <- NA
  }
  
  # Remove rows with any NA
  filtered_data <- data[complete.cases(data), ]
  final_rows <- nrow(filtered_data)
  removed_count <- initial_rows - final_rows
  
  return(list(filtered_data = filtered_data, removed_count = removed_count, outlier_info = outlier_info))
}

# Function to remove global outliers (both top and bottom percentages)
remove_global_outliers <- function(data, top_percent = 0.05, bottom_percent = 0.05, threshold_fraction = 0.5) {
  signal_cols <- 4:18
  initial_rows <- nrow(data)
  
  # Detect outliers in each column (both top and bottom)
  is_outlier <- function(column) {
    upper_threshold <- quantile(column, probs = 1 - top_percent, na.rm = TRUE)
    return(column > upper_threshold)
  }
  
  # Apply outlier detection to each signal column
  outlier_matrix <- data[, signal_cols] %>%
    apply(2, is_outlier)
  
  # Count how many samples mark each peak as an outlier
  outlier_count <- rowSums(outlier_matrix, na.rm = TRUE)
  outlier_fraction <- outlier_count / length(signal_cols)
  
  # Filter peaks where outlier fraction exceeds the threshold
  filtered_data <- data[outlier_fraction <= threshold_fraction, ]
  final_rows <- nrow(filtered_data)
  removed_count <- initial_rows - final_rows
  return(list(filtered_data = filtered_data, removed_count = removed_count))
}

# Step 1: Apply normalization
normalized_count <- normalize_tpm(count)

# Step 2: Store data for plotting before removing any outliers
all_data_no_outliers <- do.call(rbind, lapply(4:18, function(i) {
  data.frame(value = log(normalized_count[, i] + 1), sample = colnames(normalized_count)[i])
}))

# Step 3: Remove local outliers
local_outlier_result <- remove_local_outliers(normalized_count)
local_filtered_count <- local_outlier_result$filtered_data
local_removed_count <- local_outlier_result$removed_count
local_outlier_info <- local_outlier_result$outlier_info

# Save the local outlier details to a file
local_outlier_file <- "./results/local_outliers_info.txt"
write.table(local_outlier_info, file = local_outlier_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Print confirmation
cat("Local outliers information saved to:", local_outlier_file, "\n")

# Step 4: Remove global outliers
global_outlier_result <- remove_global_outliers(local_filtered_count, top_percent = 0.05, bottom_percent = 0.05, threshold_fraction = 0.5)
filtered_count <- global_outlier_result$filtered_data
global_removed_count <- global_outlier_result$removed_count

# Step 5: Store data for plotting after removing global outliers
all_data_with_outliers_removed <- do.call(rbind, lapply(4:18, function(i) {
  data.frame(value = log(filtered_count[, i] + 1), sample = colnames(filtered_count)[i])
}))

# Save normalized and filtered data (including genomic position) to a new file
output_file <- "./results/peaks_quan_tpm_filtered.txt"
write.table(filtered_count, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Save removal statistics to a separate file
stats_file <- "./results/removed_peaks_statistics.txt"
write.table(data.frame(
  Step = c("Local Outliers", "Global Outliers"),
  RemovedPeaks = c(local_removed_count, global_removed_count)
), file = stats_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Create a boxplot for all samples before removing any outliers
p1 <- ggplot(all_data_no_outliers, aes(x = sample, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12)) +
  ggtitle("Boxplots for All Samples (No Outliers Removed)") +
  ylab("Log(TPM + 1)") +
  xlab("Sample")
ggsave(file.path(output_dir, "boxplots_no_outliers_combined.pdf"), p1, width = 12, height = 6)

# Create a boxplot for all samples after removing outliers
p2 <- ggplot(all_data_with_outliers_removed, aes(x = sample, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12)) +
  ggtitle("Boxplots for All Samples (Outliers Removed)") +
  ylab("Log(TPM + 1)") +
  xlab("Sample")
ggsave(file.path(output_dir, "boxplots_outliers_removed_combined.pdf"), p2, width = 12, height = 6)

# Print success message
cat("TPM normalization, outlier removal, and visualization completed.\n")
cat("Filtered data (with position) saved to:", output_file, "\n")
cat("Statistics saved to:", stats_file, "\n")
cat("Local outliers removed:", local_removed_count, "\n")
cat("Global outliers removed:", global_removed_count, "\n")
