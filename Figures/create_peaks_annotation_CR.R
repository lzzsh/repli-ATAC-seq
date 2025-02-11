# Set working directory
setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
count <- read.table("./peaks_re_quan_CR.txt")

# Add column names
colnames(count) <- c("chr", "start", "end",
                     "xw11-OE-ES", "xw11-OE-MS", "xw11-OE-LS",
                     "xw11-CR-ES", "xw11-CR-MS", "xw11-CR-LS",
                     "x39-CR-ES", "x39-CR-MS", "x39-CR-LS",
                     "ZH11-3-ES", "ZH11-3-MS", "ZH11-3-LS",
                     "ZH11-4-ES", "ZH11-4-MS", "ZH11-4-LS")

# Directory to save the results
output_dir <- "."
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

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

# Apply normalization
normalized_count <- normalize_tpm(count)

# Prepare data for plotting
all_data <- do.call(rbind, lapply(4:18, function(i) {
  data.frame(value = log(normalized_count[, i] + 1), sample = colnames(normalized_count)[i])
}))

# Save normalized data to a new file
output_file <- "./peaks_re_quan_tpm_normalized_CR.txt"
write.table(normalized_count[,1:18], file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Create a boxplot for all samples
p <- ggplot(all_data, aes(x = sample, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.colour = "black", outlier.size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12)) +
  ggtitle("Boxplots for All Samples (Normalized Data)") +
  ylab("Log(TPM + 1)") +
  xlab("Sample")
ggsave(file.path(output_dir, "boxplots_normalized_combined_requan.pdf"), p, width = 12, height = 6)

# Print success message
cat("TPM normalization and visualization completed.\n")
cat("Normalized data saved to:", output_file, "\n")

# Step 1: Load and preprocess data
input_data <- read.table("./peaks_re_quan_tpm_normalized_CR.txt", header = TRUE) %>%
  mutate(across(where(is.numeric), ~ . + 1e-6))  # Avoid division by zero

# Step 2: Standardize signals by ZH11.3 and ZH11.4 for each phase
batches <- c("xw11.OE", "xw11.CR", "x39.CR")
normalized_data <- input_data

# Loop through each batch to normalize signals
for (batch in batches) {
  for (phase in c("ES", "MS", "LS")) {
    col <- paste0(batch, ".", phase)
    zh11_avg <- rowMeans(input_data[, c(paste0("ZH11.3.", phase), paste0("ZH11.4.", phase))])
    norm_col <- paste0(col, "_norm")
    normalized_data[[norm_col]] <- ifelse((normalized_data[[col]] / zh11_avg) > 1,
                                          normalized_data[[col]] / zh11_avg, 0)
  }
}

# Step 3: Classify regions based on signals
classified_data <- normalized_data %>%
  rowwise() %>%
  mutate(
    xw11.OE_phase = {
      signals <- c(ES = xw11.OE.ES_norm, MS = xw11.OE.MS_norm, LS = xw11.OE.LS_norm)
      if (all(signals == 0)) "Non-replication" else paste(names(signals)[signals >= 0.9 * max(signals)], collapse = "")
    },
    xw11.CR_phase = {
      signals <- c(ES = xw11.CR.ES_norm, MS = xw11.CR.MS_norm, LS = xw11.CR.LS_norm)
      if (all(signals == 0)) "Non-replication" else paste(names(signals)[signals >= 0.9 * max(signals)], collapse = "")
    },
    x39.CR_phase = {
      signals <- c(ES = x39.CR.ES_norm, MS = x39.CR.MS_norm, LS = x39.CR.LS_norm)
      if (all(signals == 0)) "Non-replication" else paste(names(signals)[signals >= 0.9 * max(signals)], collapse = "")
    }
  ) %>%
  ungroup() %>%
  rowwise() %>%
  ungroup() 

classified_data$start <- format(classified_data$start, scientific = FALSE, trim = TRUE)
classified_data$end <- format(classified_data$end, scientific = FALSE, trim = TRUE)
write.table(classified_data, "/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT_CR.gff3", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)


# Step 4: Generate GTF-like format using pre-computed phase
classified_data$final_phase <- sapply(classified_data$final_phase, function(x) gsub("S", "", x))
peaks_reads_ZH11 <- classified_data %>%
  mutate(annotation = case_when(
    final_phase == "E" ~ "Name=ES;color=#2250F1;",
    final_phase == "EM" ~ "Name=ESMS;color=#28C5CC;",
    final_phase == "M" ~ "Name=MS;color=#1A8A12;",
    final_phase == "ML" ~ "Name=MSLS;color=#FFFD33;",
    final_phase == "L" ~ "Name=LS;color=#FB0018;",
    final_phase == "EL" ~ "Name=ESLS;color=#EA3CF2;",
    final_phase == "EML" ~ "Name=ESMSLS;color=#FAB427;",
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
peaks_reads_ZH11$start <- format(peaks_reads_ZH11$start, scientific = FALSE, trim = TRUE)
peaks_reads_ZH11$end <- format(peaks_reads_ZH11$end, scientific = FALSE, trim = TRUE)
write.table(peaks_reads_ZH11, "ZH11.gff3", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)

classified_data$start <- format(classified_data$start, scientific = FALSE, trim = TRUE)
classified_data$end <- format(classified_data$end, scientific = FALSE, trim = TRUE)
selected_col <- c("chr","start","end","final_phase")
write.table(classified_data[,selected_col], "/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT.gff3", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)

cat("GTF-like file saved as 'ZH11.gff3'.\n")

# Step 6: Save results
write.csv(classified_data, "replication_classification_results.csv", row.names = FALSE)