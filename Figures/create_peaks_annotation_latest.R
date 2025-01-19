# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/")

# Step 1: Load and preprocess data
input_data <- read.table("./peaks_re_quan_tpm_normalized.txt", header = TRUE) %>%
  select(-starts_with("ZH11.3"), -starts_with("ZH11.4")) %>%  # Remove unnecessary columns
  mutate(across(where(is.numeric), ~ . + 1e-6))  # Add small value to avoid division by zero

# Step 2: Standardize ES, MS, LS columns using "ZH11.2.G1"
standardized_data <- input_data %>%
  mutate(across(matches("\\.ES$|\\.MS$|\\.LS$"),
                list(norm = ~ . / ZH11.2.G1),
                .names = "{col}_norm"))

# Step 3: Determine replication phases for each batch
classified_data <- standardized_data %>%
  rowwise() %>%
  mutate(
    # Calculate phases for each batch
    NIP.1_phase = {
      signals <- c(ES = NIP.1.ES_norm, MS = NIP.1.MS_norm, LS = NIP.1.LS_norm)
      threshold <- 0.8 * max(signals)
      paste(names(signals)[signals >= threshold], collapse = "")
    },
    ZH11.1_phase = {
      signals <- c(ES = ZH11.1.ES_norm, MS = ZH11.1.MS_norm, LS = ZH11.1.LS_norm)
      threshold <- 0.8 * max(signals)
      paste(names(signals)[signals >= threshold], collapse = "")
    },
    ZH11.2_phase = {
      signals <- c(ES = ZH11.2.ES_norm, MS = ZH11.2.MS_norm, LS = ZH11.2.LS_norm)
      threshold <- 0.8 * max(signals)
      paste(names(signals)[signals >= threshold], collapse = "")
    }
  ) %>%
  ungroup()

# Step 4: Summarize final classification by voting
classified_data <- classified_data %>%
  rowwise() %>%
  mutate(
    final_phase = names(sort(table(c(NIP.1_phase, ZH11.1_phase, ZH11.2_phase)), decreasing = TRUE)[1])
  ) %>%
  ungroup()

# Step 5: Save results
write.csv(classified_data, "replication_classification_results.csv", row.names = FALSE)