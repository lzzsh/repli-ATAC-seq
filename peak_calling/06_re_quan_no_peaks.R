library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(limma)

# Set the working directory to the location of the result files
setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/")

# Read the data from the peaks_masked_counts.txt file
wt <- read.table("./peaks_masked_counts.txt", header = FALSE)

# Set column names for the data
colnames(wt) <- c("chr", "start", "end",
                  "ZH11-2-G1", "ZH11-2-ES", "ZH11-2-MS", "ZH11-2-LS",
                  "ZH11-1-G1", "ZH11-1-ES", "ZH11-1-MS", "ZH11-1-LS",
                  "NIP-1-G1", "NIP-1-ES", "NIP-1-MS", "NIP-1-LS",
                  "xw11-OE-ES", "xw11-OE-MS", "xw11-OE-LS",
                  "xw11-CR-ES-1", "xw11-CR-MS", "xw11-CR-LS",
                  "x39-CR-ES", "x39-CR-MS", "x39-CR-LS",
                  "ZH11-3-ES", "ZH11-3-MS", "ZH11-3-LS",
                  "ZH11-4-ES", "ZH11-4-MS", "ZH11-4-LS",
                  "xw11-CR-ES-2")
# Add RT and cuttag columns (assuming these are present in another dataset called 'data')
wt$RT <- data$RT
wt$cuttag <- data$cuttag

# Replace all negative values in the data with 0
wt[wt < 0] <- 0 

# Select columns corresponding to signal values (excluding non-numeric columns)
signal_cols <- 4:(ncol(wt) - 2)  # Assume the last two columns are RT and cuttag

# Apply log2 transformation with +1 to the signal values
wt[signal_cols] <- log2(wt[signal_cols] + 1)

# Perform quantile normalization across arrays (samples)
wt[signal_cols] <- normalizeBetweenArrays(as.matrix(wt[signal_cols]), method = "quantile")

# Convert the 'cuttag' column to a factor
wt$cuttag <- as.factor(wt$cuttag)

# Step to filter "overlap" points from the data and exclude non-replication and unknown RT values
overlap_points <- subset(wt, cuttag == "overlap" & !RT %in% c("Non-replication", "unknown"))

# Calculate the number of points above and below the y=x reference line
above_yx <- sum(overlap_points$`xw11-CR-ES-2` > overlap_points$`ZH11-2-ES`)  # Points above the line (y > x)
below_yx <- sum(overlap_points$`xw11-CR-ES-2` < overlap_points$`ZH11-2-ES`)  # Points below the line (y < x)

# Calculate the total number of overlap points
total_overlap <- nrow(overlap_points)

# Calculate the percentage of points above and below the y=x reference line
above_ratio <- above_yx / total_overlap * 100
below_ratio <- below_yx / total_overlap * 100

# Create a scatter plot to visualize the data
p <- ggplot(wt) +
  # Plot gray points for the entire data
  geom_point(aes(x = `ZH11-2-ES`, y = `xw11-CR-ES-2`), color = "grey", alpha = 0.5, size = 1.5) +
  # Plot red points for the overlap points
  geom_point(data = overlap_points, aes(x = `ZH11-2-ES`, y = `xw11-CR-ES-2`), 
             color = "red", size = 1.8, alpha = 0.8) +
  # Add a dashed reference line for y = x
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
  # Set aspect ratio of x and y axes to 1:1 for consistent scaling
  coord_fixed(ratio = 1) +
  labs(
    x = "ZH11-2-ES (log2 quantile norm)",
    y = "xw11-CR-ES (log2 quantile norm)",
    title = "Scatter Plot with y=x Line and Overlap Analysis"
  ) +
  theme_minimal()

# Print the overlap ratio information
cat(sprintf("Red points: %.2f%% are above the y=x line, %.2f%% are below the y=x line.\n", above_ratio, below_ratio))

# Display the plot
print(p)
