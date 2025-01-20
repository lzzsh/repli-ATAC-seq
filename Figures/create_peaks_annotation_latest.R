# Load necessary libraries
library(dplyr)
library(tidyr)

setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/")

# Step 1: Load and preprocess data
input_data <- read.table("./peaks_re_quan_tpm_normalized.txt", header = TRUE) %>%
  select(-starts_with("ZH11.3"), -starts_with("ZH11.4")) %>%
  mutate(across(where(is.numeric), ~ . + 1e-6))  # Avoid division by zero

# Step 2: Sum up G1, ES, MS, and LS signals across all batches
summed_data <- input_data %>%
  mutate(
    total_G1 = rowSums(across(ends_with(".G1"))),
    total_ES = rowSums(across(ends_with(".ES"))),
    total_MS = rowSums(across(ends_with(".MS"))),
    total_LS = rowSums(across(ends_with(".LS")))
  )

# Step 3: Normalize signals and apply thresholding (<1 becomes 0)
classified_data <- summed_data %>%
  mutate(
    norm_ES = ifelse(total_ES / total_G1 >= 1, total_ES / total_G1, 0),
    norm_MS = ifelse(total_MS / total_G1 >= 1, total_MS / total_G1, 0),
    norm_LS = ifelse(total_LS / total_G1 >= 1, total_LS / total_G1, 0)
  ) %>%
  rowwise() %>%
  mutate(
    phase = {
      signals <- c(ES = norm_ES, MS = norm_MS, LS = norm_LS)
      
      if (all(signals == 0)) {
        "Non-replication"
      } else {
        valid_phases <- names(signals)[signals >= 0.9 * max(signals)]
        paste(valid_phases, collapse = "")
      }
    }
  ) %>%
  ungroup() %>%
  filter(phase != "Non-replication")  # Remove Non-replication rows

# Step 4: Generate GTF-like format using pre-computed phase
classified_data$phase <- sapply(classified_data$phase, function(x) gsub("S", "", x))
peaks_reads_ZH11 <- classified_data %>%
  mutate(annotation = case_when(
    phase == "E" ~ "Name=ES;color=#2250F1;",
    phase == "EM" ~ "Name=ESMS;color=#28C5CC;",
    phase == "M" ~ "Name=MS;color=#1A8A12;",
    phase == "ML" ~ "Name=MSLS;color=#FFFD33;",
    phase == "L" ~ "Name=LS;color=#FB0018;",
    phase == "EL" ~ "Name=ESLS;color=#EA3CF2;",
    phase == "EML" ~ "Name=ESMSLS;color=#FAB427;",
    TRUE ~ ""
  )) %>%
  mutate(
    source = ".",
    score = ".",
    strand = ".",
    frame = ".",
    feature = "peaks"
  ) %>%
  select(chr, source, feature, start, end, score, strand, frame, annotation)

# Step 5: Save GTF file
write.table(peaks_reads_ZH11, "ZH11.gff3", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)

cat("GTF-like file saved as 'ZH11.gff3'.\n")


# # Load necessary libraries
# library(dplyr)
# library(tidyr)
# 
# setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/")
# 
# # Step 1: Load and preprocess data
# input_data <- read.table("./peaks_re_quan_tpm_normalized.txt", header = TRUE) %>%
#   select(-starts_with("ZH11.3"), -starts_with("ZH11.4")) %>%
#   mutate(across(where(is.numeric), ~ . + 1e-6))  # Avoid division by zero
# 
# # Step 2: Standardize signals by G1 for each batch
# batches <- c("NIP.1", "ZH11.1", "ZH11.2")
# normalized_data <- input_data
# 
# # Loop through each batch to normalize signals
# for (batch in batches) {
#   g1_col <- paste0(batch, ".G1")
#   for (phase in c("ES", "MS", "LS")) {
#     col <- paste0(batch, ".", phase)
#     norm_col <- paste0(col, "_norm")
#     normalized_data[[norm_col]] <- ifelse((normalized_data[[col]] / normalized_data[[g1_col]]) > 1,
#                                           normalized_data[[col]] / normalized_data[[g1_col]], 0)
#   }
# }
# 
# # Step 3: Classify regions based on signals
# classified_data <- normalized_data %>%
#   rowwise() %>%
#   mutate(
#     NIP.1_phase = {
#       signals <- c(ES = NIP.1.ES_norm, MS = NIP.1.MS_norm, LS = NIP.1.LS_norm)
#       if (all(signals == 0)) "Non-replication" else paste(names(signals)[signals >= 0.8 * max(signals)], collapse = "")
#     },
#     ZH11.1_phase = {
#       signals <- c(ES = ZH11.1.ES_norm, MS = ZH11.1.MS_norm, LS = ZH11.1.LS_norm)
#       if (all(signals == 0)) "Non-replication" else paste(names(signals)[signals >= 0.8 * max(signals)], collapse = "")
#     },
#     ZH11.2_phase = {
#       signals <- c(ES = ZH11.2.ES_norm, MS = ZH11.2.MS_norm, LS = ZH11.2.LS_norm)
#       if (all(signals == 0)) "Non-replication" else paste(names(signals)[signals >= 0.8 * max(signals)], collapse = "")
#     }
#   ) %>%
#   ungroup() %>%
#   rowwise() %>%
#   mutate(
#     final_phase = ifelse(all(c(NIP.1_phase, ZH11.1_phase, ZH11.2_phase) == "Non-replication"),
#                          "Non-replication",
#                          names(sort(table(c(NIP.1_phase, ZH11.1_phase, ZH11.2_phase)), decreasing = TRUE)[1]))
#   ) %>%
#   ungroup()
# 
# # Step 4: Save results
# write.csv(classified_data, "replication_classification_results.csv", row.names = FALSE)


