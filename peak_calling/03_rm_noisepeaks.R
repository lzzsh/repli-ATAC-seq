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
                     "ZH11-1-ES", "ZH11-1-MS", "ZH11-1-LS",
                     "NIP-1-ES", "NIP-1-MS", "NIP-1-LS")

# Directory to save the results
output_dir <- "./results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Initialize data frames for plotting
all_data_no_outliers <- data.frame()
all_data_with_outliers_removed <- data.frame()

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

# Function to remove local outliers (within each sample) using IQR
remove_local_outliers <- function(data) {
  signal_cols <- 4:19
  initial_rows <- nrow(data)
  
  # Apply IQR-based filtering to each signal column
  data[, signal_cols] <- lapply(data[, signal_cols], function(col) {
    Q1 <- quantile(col, 0.25, na.rm = TRUE)
    Q3 <- quantile(col, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- max(0, Q1 - 1.5 * IQR)  # Ensure lower bound is not negative
    upper_bound <- Q3 + 1.5 * IQR
    col[col < lower_bound | col > upper_bound] <- NA  # Set outliers to NA
    return(col)
  })
  
  # Remove rows with any NA
  filtered_data <- data[complete.cases(data), ]
  final_rows <- nrow(filtered_data)
  removed_count <- initial_rows - final_rows
  return(list(filtered_data = filtered_data, removed_count = removed_count))
}

# Function to remove global outliers (both top and bottom percentages)
remove_global_outliers <- function(data, top_percent = 0.05, bottom_percent = 0.05, threshold_fraction = 0.5) {
  signal_cols <- 4:19
  initial_rows <- nrow(data)
  
  # Detect outliers in each column (both top and bottom)
  is_outlier <- function(column) {
    upper_threshold <- quantile(column, probs = 1 - top_percent, na.rm = TRUE)
    lower_threshold <- quantile(column, probs = bottom_percent, na.rm = TRUE)
    return(column > upper_threshold | column < lower_threshold)
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
all_data_no_outliers <- do.call(rbind, lapply(4:19, function(i) {
  data.frame(value = log(normalized_count[, i] + 1), sample = colnames(normalized_count)[i])
}))

# Step 3: Remove local outliers
local_outlier_result <- remove_local_outliers(normalized_count)
local_filtered_count <- local_outlier_result$filtered_data
local_removed_count <- local_outlier_result$removed_count

# Step 4: Remove global outliers
global_outlier_result <- remove_global_outliers(local_filtered_count, top_percent = 0.05, bottom_percent = 0.05, threshold_fraction = 0.5)
filtered_count <- global_outlier_result$filtered_data
global_removed_count <- global_outlier_result$removed_count

# Step 5: Store data for plotting after removing global outliers
all_data_with_outliers_removed <- do.call(rbind, lapply(4:19, function(i) {
  data.frame(value = log(filtered_count[, i] + 1), sample = colnames(filtered_count)[i])
}))

# Save normalized and filtered data (including genomic position) to a new file
output_file <- "./results/peaks_quan_tpm_filtered.txt"
write.table(filtered_count, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

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