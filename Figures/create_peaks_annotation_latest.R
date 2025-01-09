# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Load Input Data
input_data <- read.csv("input_table.csv")  # Input file with columns: "chr", "peak", "ES_rep1", ..., "G1_rep1"

# Step 2: Calculate Average Signal for Each Replication Phase
averaged_data <- input_data %>%
  rowwise() %>%
  mutate(
    ES_mean = mean(c_across(starts_with("ES_rep")), na.rm = TRUE),
    MS_mean = mean(c_across(starts_with("MS_rep")), na.rm = TRUE),
    LS_mean = mean(c_across(starts_with("LS_rep")), na.rm = TRUE),
    G1_mean = mean(c_across(starts_with("G1_rep")), na.rm = TRUE)
  ) %>%
  ungroup()

# Step 3: Normalize Signal (Rcwt)
normalized_data <- averaged_data %>%
  mutate(
    ES_normalized = ES_mean / G1_mean,
    MS_normalized = MS_mean / G1_mean,
    LS_normalized = LS_mean / G1_mean
  )

# Step 4: Determine Replication Thresholds (mTc and ˆTc) for Each Chromosome
threshold_results <- normalized_data %>%
  group_by(chr) %>%
  summarise(
    mTc = {
      replication_signals <- c(ES_normalized, MS_normalized, LS_normalized)
      sorted_signals <- sort(replication_signals, decreasing = FALSE)
      coverage_rate <- seq_along(sorted_signals) / length(sorted_signals)
      coverage_diff <- c(0, diff(coverage_rate) / diff(sorted_signals))
      sorted_signals[which.max(abs(coverage_diff))]
    },
    hat_Tc = {
      replication_signals <- c(ES_normalized, MS_normalized, LS_normalized)
      sorted_signals <- sort(replication_signals, decreasing = FALSE)
      coverage_rate <- seq_along(sorted_signals) / length(sorted_signals)
      coverage_diff <- c(0, diff(coverage_rate) / diff(sorted_signals))
      valid_indices <- which(sorted_signals < mTc & abs(coverage_diff) < 0.1)
      if (length(valid_indices) > 0) sorted_signals[max(valid_indices)] else mTc
    },
    .groups = "drop"
  )

# Step 5: Apply Replication Thresholds to Classify Replication Signal
classified_data <- normalized_data %>%
  left_join(threshold_results, by = "chr") %>%
  mutate(
    ES_classification = ifelse(ES_normalized > hat_Tc, "Replicating", "Non-replicating"),
    MS_classification = ifelse(MS_normalized > hat_Tc, "Replicating", "Non-replicating"),
    LS_classification = ifelse(LS_normalized > hat_Tc, "Replicating", "Non-replicating")
  )

# Step 6: Assign Multiple Replication Phases
# Determine all replication phases where the signal exceeds the threshold
classified_data <- classified_data %>%
  rowwise() %>%
  mutate(
    replication_phases = paste(
      c(
        if (ES_normalized > hat_Tc) "Early S-phase" else NULL,
        if (MS_normalized > hat_Tc) "Middle S-phase" else NULL,
        if (LS_normalized > hat_Tc) "Late S-phase" else NULL
      ),
      collapse = ", "
    ),
    replication_phases = ifelse(replication_phases == "", "Non-replicating", replication_phases)
  ) %>%
  ungroup()

# Step 7: Visualize Classification Results
plot_signal_distribution <- function(data, phase, title) {
  ggplot(data, aes_string(x = paste0(phase, "_normalized"), fill = paste0(phase, "_classification"))) +
    geom_histogram(binwidth = 0.1, alpha = 0.8, position = "stack") +
    labs(title = title, x = paste0("Normalized Signal (", phase, ")"), y = "Count") +
    theme_minimal()
}

plot_signal_distribution(classified_data, "ES", "Early S-phase Signal Distribution")
plot_signal_distribution(classified_data, "MS", "Middle S-phase Signal Distribution")
plot_signal_distribution(classified_data, "LS", "Late S-phase Signal Distribution")

# Step 8: Save Results
write.csv(classified_data, "replication_classification_multiple_phases.csv", row.names = FALSE)