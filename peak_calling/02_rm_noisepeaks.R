# Set working directory
setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# List all the .txt files in the directory
files <- list.files("./", pattern = "*.txt", full.names = TRUE)

# Extract batch (before the last '_') and period (after the last '_') for sorting
file_info <- data.frame(
  file = files,
  batch = sapply(basename(files), function(x) sub("_(ES|MS|LS).*", "", x)),  # Extract everything before the last '_period'
  period = sapply(basename(files), function(x) sub(".*_(ES|MS|LS)\\..*", "\\1", x))  # Extract the text after the last '_'
)

# Define the order for periods (ES, MS, LS)
period_order <- c("ES", "MS", "LS")
file_info$period <- factor(file_info$period, levels = period_order)

# Sort files by batch and then period
file_info <- file_info %>% arrange(batch, period)
files <- file_info$file

# Update file names for labeling
file_info$label <- basename(file_info$file)

# Create empty data frames to store the data for ggplot
all_data_no_outliers <- data.frame()
all_data_with_outliers_removed <- data.frame()

# Directory to save the files with outliers removed
output_dir <- "./results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Data frame to store the outlier removal statistics
outlier_stats <- data.frame(Sample = character(), Total = numeric(), Removed = numeric(), Percentage = numeric())

# Phase: Using boxplot method to remove outliers
for (i in seq_along(files)) {
  file <- files[i]
  label <- basename(file)  # Use basename for consistent labeling
  
  # Read the data
  a <- read.table(file, sep = "\t", header = FALSE)  # Adjust `header` if files have column names
  
  # Log transformation (column 11)
  a[, 11] <- a[, 11] + 1
  data <- log(a[, 11]) 
  
  # Store the data for boxplot (without removing outliers)
  all_data_no_outliers <- rbind(all_data_no_outliers, 
                                data.frame(value = data, sample = label))
  
  # Remove outliers using IQR (boxplot method)
  Q1 <- quantile(data, 0.25)  # First quartile
  Q3 <- quantile(data, 0.75)  # Third quartile
  IQR <- Q3 - Q1             # Interquartile range
  
  # Define bounds for outliers
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # Filter data within bounds
  data_trimmed <- a[data >= lower_bound & data <= upper_bound, ]
  
  # Store the trimmed data for boxplot after removing outliers
  all_data_with_outliers_removed <- rbind(all_data_with_outliers_removed, 
                                          data.frame(value = log(data_trimmed[, 11]), sample = label))
  
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
  
  # Save the trimmed data into a new file
  output_file <- file.path(output_dir, paste0("filtered_", label))
  write.table(data_trimmed, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Sort the outlier stats by batch and period
outlier_stats <- outlier_stats %>%
  mutate(
    period = factor(sub(".*_(.*)\\..*", "\\1", Sample), levels = period_order),  # Set the order for periods
    batch = factor(sub("_(ES|MS|LS).*", "", Sample), levels = unique(sub("_(ES|MS|LS).*", "", Sample)))  # Convert batch to factor
  ) %>%
  arrange(batch, period)  # Sort by batch and period

# Print the sorted outlier removal stats
print(outlier_stats)

# Save the sorted outlier stats to a CSV file
write.csv(outlier_stats, file = file.path(output_dir, "outlier_removal_stats_sorted.csv"), row.names = FALSE)

# Create a boxplot without removing outliers
p1 <- ggplot(all_data_no_outliers, aes(x = factor(sample, levels = file_info$label), y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12)) +
  ggtitle("Boxplots (No Outliers Removed)") +
  ylab("Log Transformed Data") +
  xlab("Sample")
ggsave(file.path(output_dir, "boxplots_no_outliers.pdf"), p1, width = 10, height = 6)

# Create a boxplot after removing outliers
p2 <- ggplot(all_data_with_outliers_removed, aes(x = factor(sample, levels = file_info$label), y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12)) +
  ggtitle("Boxplots (Outliers Removed)") +
  ylab("Log Transformed Data") +
  xlab("Sample")
ggsave(file.path(output_dir, "boxplots_outliers_removed.pdf"), p2, width = 10, height = 6)

# Create a density plot for all files (no outliers removed, sorted by files)
p3 <- ggplot(all_data_no_outliers, aes(x = value, fill = sample)) +
  geom_density(alpha = 0.5) +
  ggtitle("Density Plots (No Outliers Removed)") +
  xlab("Log Transformed Data") +
  ylab("Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12)) +
  facet_wrap(~factor(sample, levels = file_info$label), scales = "free", ncol = 5)  # Facet by sample
ggsave(file.path(output_dir, "density_plots_no_outliers.pdf"), p3, width = 12, height = 8)