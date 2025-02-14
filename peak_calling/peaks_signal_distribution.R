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
  
  # Data frame to store IQR bounds for each sample
  iqr_bounds <- data.frame(Sample = character(), Q1 = numeric(), Q3 = numeric(),
                           LowerBound = numeric(), UpperBound = numeric(), stringsAsFactors = FALSE)
  
  # Apply IQR-based filtering to each signal column
  for (col_idx in signal_cols) {
    col_name <- colnames(data)[col_idx]
    Q1 <- quantile(data[[col_idx]], 0.25, na.rm = TRUE)
    Q3 <- quantile(data[[col_idx]], 0.75, na.rm = TRUE)
    IQR_value <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR_value
    upper_bound <- Q3 + 1.5 * IQR_value
    
    # Store IQR bounds
    iqr_bounds <- rbind(iqr_bounds, data.frame(Sample = col_name, Q1 = Q1, Q3 = Q3,
                                               LowerBound = lower_bound, UpperBound = upper_bound))
    
    # Print IQR bounds for each sample
    cat("Sample:", col_name, "\n",
        "Q1:", Q1, "| Q3:", Q3, "\n",
        "Lower Bound:", lower_bound, "| Upper Bound:", upper_bound, "\n\n")
    
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
  
  # Save IQR bounds to a file
  write.table(iqr_bounds, file = "./results/IQR_bounds.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save outlier information to a file
  write.table(outlier_info, file = "./results/local_outliers_info.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("IQR bounds saved to: ./results/IQR_bounds.txt\n")
  cat("Local outliers information saved to: ./results/local_outliers_info.txt\n")
  
  # Remove rows with any NA
  filtered_data <- data[complete.cases(data), ]
  final_rows <- nrow(filtered_data)
  removed_count <- initial_rows - final_rows
  
  return(list(filtered_data = filtered_data, removed_count = removed_count, 
              outlier_info = outlier_info, iqr_bounds = iqr_bounds))
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

# signal histplot
org <- read.table("./organelle_peaks.bed")
org <- org %>% unique()
colnames(org)[1:3] <- c("chr","start","end")
org_no_outliers <- read.table("./organelle_peaks_no_outliers.bed")
org_no_outliers <- org_no_outliers %>% unique()
colnames(org_no_outliers)[c(1,4,5)] <- c("chr","start","end")

# Select the first three columns from org and org_no_outliers
org_subset <- org %>% select(1:3)
org_no_outliers_subset <- org_no_outliers %>% select(1,4,5)

# Filter normalized_count to keep rows that match org's first three columns
filtered_by_org <- normalized_count %>%
  semi_join(org_subset, by = colnames(org_subset))

# Filter normalized_count to keep rows that match org_no_outliers' first three columns
filtered_by_org_no_outliers <- normalized_count %>%
  semi_join(org_no_outliers_subset, by = colnames(org_no_outliers_subset))

col_names = c(
             "ZH11-3-ES", "ZH11-3-MS", "ZH11-3-LS",
             "ZH11-4-ES", "ZH11-4-MS", "ZH11-4-LS",
             "ZH11-2-ES", "ZH11-2-MS", "ZH11-2-LS",
             "ZH11-1-ES", "ZH11-1-MS", "ZH11-1-LS",
             "NIP-1-ES", "NIP-1-MS", "NIP-1-LS")

count <- filtered_by_org

# Select signal columns (assuming the first three columns are chr, start, and end)
signal_cols <- 4:(ncol(count)-1)

# Convert the data to long format (tidy data)
long_data <- count %>%
  select(all_of(signal_cols)) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Signal")

# Keep only signal values greater than 50 and apply log2 transformation
long_data <- long_data %>%
  mutate(Log2_Signal = log2(Signal))

# Plot histogram of Log2(Signal) across all samples
ggplot(long_data, aes(x = Log2_Signal, fill = Sample, color = Sample)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  labs(title = "Histogram of Log2(Signal) Across All Samples",
       x = "Log2(Signal)", 
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

# Plot density plot of Log2(Signal) across all samples
ggplot(long_data, aes(x = Log2_Signal, fill = Sample, color = Sample)) +
  geom_density(alpha = 0.3) +  # Draw density curves with transparency
  labs(title = "Density Plot of Log2(Signal) Across All Samples",
       x = "Log2(Signal)", 
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title