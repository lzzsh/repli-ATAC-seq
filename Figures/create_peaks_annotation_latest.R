# Load necessary libraries
library(dplyr)
library(ggplot2)
library(scales)

# Step 1: Load Biological Replicates Data
# Read biological replicates for S-phase
s_rep1 <- read.table("peak_coverage_S_rep1.bed", header = FALSE, col.names = c("chr", "start", "end", "peak_id", "coverage_rep1"))
s_rep2 <- read.table("peak_coverage_S_rep2.bed", header = FALSE, col.names = c("chr", "start", "end", "peak_id", "coverage_rep2"))

# Read biological replicates for G1-phase
g1_rep1 <- read.table("peak_coverage_G1_rep1.bed", header = FALSE, col.names = c("chr", "start", "end", "peak_id", "coverage_rep1"))
g1_rep2 <- read.table("peak_coverage_G1_rep2.bed", header = FALSE, col.names = c("chr", "start", "end", "peak_id", "coverage_rep2"))

# Step 2: Merge Replicates and Calculate Mean Coverage
# Merge S-phase replicates
s_combined <- s_rep1 %>%
  inner_join(s_rep2, by = c("chr", "start", "end", "peak_id")) %>%
  mutate(coverage_mean_S = rowMeans(select(., coverage_rep1, coverage_rep2)))

# Merge G1-phase replicates
g1_combined <- g1_rep1 %>%
  inner_join(g1_rep2, by = c("chr", "start", "end", "peak_id")) %>%
  mutate(coverage_mean_G1 = rowMeans(select(., coverage_rep1, coverage_rep2)))

# Merge S-phase and G1-phase mean coverage data
merged_data <- s_combined %>%
  select(chr, start, end, peak_id, coverage_mean_S) %>%
  inner_join(g1_combined %>% select(chr, start, end, peak_id, coverage_mean_G1), by = c("chr", "start", "end", "peak_id"))

# Step 3: Calculate Normalized Signal (Signal Ratio)
merged_data <- merged_data %>%
  mutate(signal_ratio = coverage_mean_S / coverage_mean_G1)

# View the first few rows of normalized data
head(merged_data)

# Step 4: Dynamically Determine Optimal Signal Threshold
# Sort data by signal ratio
sorted_data <- merged_data %>%
  arrange(signal_ratio)

# Calculate coverage rate and rate of change in coverage
thresholds <- sorted_data %>%
  mutate(coverage_rate = row_number() / n()) %>%
  mutate(coverage_diff = c(0, diff(coverage_rate) / diff(signal_ratio)))

# Identify the optimal threshold where the rate of change is maximum
optimal_threshold <- thresholds %>%
  filter(coverage_diff == max(coverage_diff, na.rm = TRUE)) %>%
  pull(signal_ratio)

# Print the optimal signal threshold
cat("Optimal Signal Threshold:", optimal_threshold, "\n")

# Step 5: Classify Replication Phases
# Set thresholds for classification
early_threshold <- optimal_threshold
middle_threshold <- 2 * optimal_threshold

# Classify each peak into replication phases
merged_data <- merged_data %>%
  mutate(replication_phase = case_when(
    signal_ratio >= middle_threshold ~ "Late",
    signal_ratio >= early_threshold & signal_ratio < middle_threshold ~ "Middle",
    signal_ratio < early_threshold ~ "Early",
    TRUE ~ "Non-replicating"
  ))

# View classified data
head(merged_data)

# Step 6: Visualize Signal Distribution and Replication Phases
# Plot histogram of normalized signal ratio with replication phase classification
ggplot(merged_data, aes(x = signal_ratio, fill = replication_phase)) +
  geom_histogram(binwidth = 0.1, alpha = 0.8, position = "stack") +
  scale_x_log10() +
  labs(title = "Replication Phase Distribution",
       x = "Normalized Signal Ratio",
       y = "Peak Count") +
  theme_minimal()

# Step 7: Save Results
# Save classified data to a file
write.table(merged_data, "replication_phase_classification_with_replicates.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Final classified file format:
# chr    start    end    peak_id    coverage_mean_S    coverage_mean_G1    signal_ratio    replication_phase
# chr1   1000     2000   peak_1    150                100                 1.5             Middle
# chr1   2000     3000   peak_2    200                50                  4.0             Late