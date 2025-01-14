# Set working directory
setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# List all the .txt files in the directory
files <- list.files("./", pattern = "*.txt", full.names = TRUE)

# Extract batch and period information for sorting
file_info <- data.frame(
  file = files,
  batch = sapply(basename(files), function(x) sub("_(ES|MS|LS).*", "", x)),  # Extract batch
  period = sapply(basename(files), function(x) sub(".*_(ES|MS|LS)\\..*", "\\1", x))  # Extract period
)

# Define the order for periods
period_order <- c("ES", "MS", "LS")
file_info$period <- factor(file_info$period, levels = period_order)

# Sort files by batch and period
file_info <- file_info %>% arrange(batch, period)
files <- file_info$file

# Add labels to file_info for consistent x-axis
file_info$label <- sub("\\.txt$", "", basename(file_info$file))

# Directory to save the results
output_dir <- "./results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Data frame to store outlier removal statistics
outlier_stats <- data.frame(Sample = character(), Total = numeric(), Removed = numeric(), Percentage = numeric())

# Function to calculate region length
calculate_region_length <- function(data) {
  # Assuming the first three columns are: Chromosome, Start, End
  data <- data %>%
    mutate(RegionLength = V3 - V2 + 1)  # V3 = End, V2 = Start
  return(data)
}

# Function to normalize using TPM
calculate_tpm <- function(data) {
  # Calculate Reads Per Kilobase (RPK)
  data <- data %>%
    mutate(RPK = V11 / (RegionLength / 1000))  
  
  # Calculate the scaling factor
  scaling_factor <- sum(data$RPK) / 1e6
  
  # Calculate TPM
  data <- data %>%
    mutate(TPM = RPK / scaling_factor)
  
  return(data)
}

# Create empty data frames to store the data for ggplot
all_data_no_outliers <- data.frame()
all_data_with_outliers_removed <- data.frame()

# Process each file: Normalize using TPM and remove outliers
for (i in seq_along(files)) {
  file <- files[i]
  label <- sub("\\.txt$", "", basename(file))  # Extract sample name
  
  # Read the data
  a <- read.table(file, sep = "\t", header = FALSE)  # Adjust `header` if files have column names
  
  # Add region length
  a <- calculate_region_length(a)
  
  # Normalize using TPM
  a <- calculate_tpm(a)
  
  # Log transformation (TPM values)
  data <- log(a$TPM + 1)
  
  # Store the data for boxplot (without removing outliers)
  all_data_no_outliers <- rbind(all_data_no_outliers, 
                                data.frame(value = data, sample = label))  # Ensure label is stored as sample name
  
  # Remove outliers using IQR (boxplot method)
  Q1 <- quantile(data, 0.25)  # First quartile
  Q3 <- quantile(data, 0.75)  # Third quartile
  IQR <- Q3 - Q1             # Interquartile range
  
  # Define bounds for outliers
  lower_bound <- Q1 - 1.0 * IQR
  upper_bound <- Q3 + 1.0 * IQR
  
  # Filter data within bounds
  data_trimmed <- a[data >= lower_bound & data <= upper_bound, ]
  
  # Store the trimmed data for boxplot after removing outliers
  all_data_with_outliers_removed <- rbind(all_data_with_outliers_removed, 
                                          data.frame(value = log(data_trimmed$TPM + 1), sample = label))
  
  # Calculate outlier removal stats
  total_values <- nrow(a)
  removed_values <- total_values - nrow(data_trimmed)
  percentage_removed <- (removed_values / total_values) * 100
  
  # Append to the stats data frame
  outlier_stats <- rbind(outlier_stats, 
                         data.frame(Sample = label, 
                                    Total = total_values, 
                                    Removed = removed_values, 
                                    Percentage = percentage_removed))
  
  # Save the normalized and trimmed data into a new file
  output_file <- file.path(output_dir, paste0("filtered_tpm_", label, ".txt"))
  write.table(data_trimmed, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# Sort the outlier stats by batch and period
outlier_stats <- outlier_stats %>%
  mutate(
    period = factor(sub(".*_(.*)\\..*", "\\1", Sample), levels = period_order),  # Set the order for periods
    batch = factor(sub("_(ES|MS|LS).*", "", Sample), levels = unique(sub("_(ES|MS|LS).*", "", Sample)))  # Convert batch to factor
  ) %>%
  arrange(batch, period)  # Sort by batch and period

# Save the sorted outlier stats to a CSV file
write.csv(outlier_stats, file = file.path(output_dir, "outlier_removal_stats_sorted.csv"), row.names = FALSE)

# Create a combined boxplot for all samples in one figure (x-axis as sample)
p1 <- ggplot(all_data_no_outliers, aes(x = factor(sample, levels = file_info$label), y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12)) +
  ggtitle("Boxplots for All Samples (No Outliers Removed)") +
  ylab("Log(TPM + 1)") +
  xlab("Sample")
ggsave(file.path(output_dir, "boxplots_no_outliers_combined.pdf"), p1, width = 12, height = 6)

# Create a combined boxplot for all samples in one figure (x-axis as sample, outliers removed)
p2 <- ggplot(all_data_with_outliers_removed, aes(x = factor(sample, levels = file_info$label), y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12)) +
  ggtitle("Boxplots for All Samples (Outliers Removed)") +
  ylab("Log(TPM + 1)") +
  xlab("Sample")
ggsave(file.path(output_dir, "boxplots_outliers_removed_combined.pdf"), p2, width = 12, height = 6)

# Print message when processing is done
cat("TPM normalization, outlier removal, and visualization completed. Results saved in", output_dir, "\n")