setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# List all the .txt files in the directory
files <- list.files("./", pattern = "*.txt", full.names = TRUE)

# Create empty data frames to store the data for ggplot
all_data_no_outliers <- data.frame()
all_data_with_outliers_removed <- data.frame()

# Directory to save the files with outliers removed
output_dir <- "./output_no_outliers/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Phase 1: Without removing outliers (Density and Boxplot)
# Loop through each file and process
for (file in files) {
  # Read the data
  a <- read.table(file, sep = "\t", header = FALSE)  # Adjust `header` if files have column names
  
  # Log transformation (column 11)
  a[, 11] <- a[, 11] + 1
  data <- log(a[, 11]) 
  
  # Store the data for boxplot (without removing outliers)
  all_data_no_outliers <- rbind(all_data_no_outliers, 
                                data.frame(value = data, sample = basename(file)))
  
  # Remove outliers (bottom and top 0.5%)
  lower_bound <- quantile(data, 0.005)
  upper_bound <- quantile(data, 0.995)
  data_trimmed <- a[data >= lower_bound & data <= upper_bound, ]  # Filter original rows
  
  # Store the trimmed data for boxplot after removing outliers
  all_data_with_outliers_removed <- rbind(all_data_with_outliers_removed, 
                                          data.frame(value = log(data_trimmed[, 11]), sample = basename(file)))
  
  # Save the trimmed data into a new file
  output_file <- file.path(output_dir, paste0("filtered_", basename(file)))
  write.table(data_trimmed, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Create a boxplot without removing outliers (outlier size smaller, color black)
ggplot(all_data_no_outliers, aes(x = sample, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +  # Black outliers with smaller size
  theme_minimal() +  # Use minimal theme (white background)
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Rotate x-axis labels by 45 degrees
        axis.title = element_text(size = 10),  # Adjust axis title size
        plot.title = element_text(size = 12)) +  # Adjust plot title size
  ggtitle("Boxplots (No Outliers Removed)") +
  ylab("Log Transformed Data") +
  xlab("Sample")

# Create a boxplot after removing outliers (outlier size smaller, color black)
ggplot(all_data_with_outliers_removed, aes(x = sample, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +  # Black outliers with smaller size
  theme_minimal() +  # Use minimal theme (white background)
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Rotate x-axis labels by 45 degrees
        axis.title = element_text(size = 10),  # Adjust axis title size
        plot.title = element_text(size = 12)) +  # Adjust plot title size
  ggtitle("Boxplots (Outliers Removed)") +
  ylab("Log Transformed Data") +
  xlab("Sample")

# Density plots for all files (without removing outliers) in a single plot
ggplot(all_data_no_outliers, aes(x = value, fill = sample)) +
  geom_density(alpha = 0.5) +  # Semi-transparent density plots
  ggtitle("Density Plots (No Outliers Removed)") +
  xlab("Log Transformed Data") +
  ylab("Density") +
  theme_minimal() +  # Use minimal theme (white background)
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Rotate x-axis labels
        axis.title = element_text(size = 10),  # Adjust axis title size
        plot.title = element_text(size = 12)) +  # Adjust plot title size
  facet_wrap(~sample, scales = "free", ncol = 5)  # Facet by sample with separate plots but within one frame
