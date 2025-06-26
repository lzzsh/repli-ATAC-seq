# Set working directory
setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation")

# Load necessary libraries
library(dplyr)

# Step 0: Load raw counts and apply TPM normalization first
count_matrix <- read.table("repli_peaks_quan.txt")
colnames(count_matrix) <- c(
  "chr","start","end",
  "WT-1-G1", "WT-2-G1", "WT-1-ES", "WT-2-ES", "WT-1-MS", "WT-1-LS",
  "sol1_5-1-G1", "sol1_5-2-G1", "sol1_5-1-ES", "sol1_5-2-ES", "sol1_5-1-MS", "sol1_5-1-LS",
  "sol1_8-1-G1", "sol1_8-2-G1", "sol1_8-1-ES", "sol1_8-2-ES", "sol1_8-1-MS", "sol1_8-1-LS",
  "tcx2_1-1-G1", "tcx2_1-2-G1", "tcx2_1-1-ES", "tcx2_1-2-ES", "tcx2_1-1-MS", "tcx2_1-1-LS",
  "tcx2_3-1-G1", "tcx2_3-1-ES", "tcx2_3-1-MS", "tcx2_3-1-LS"
)

# TPM normalization function
normalize_tpm <- function(data) {
  data$RegionLength <- data$end - data$start + 1
  data[, 4:31] <- lapply(data[, 4:31], function(col) {
    rpk <- col / (data$RegionLength / 1000)
    scaling_factor <- sum(rpk) / 1e6
    tpm <- rpk / scaling_factor
    return(tpm)
  })
  return(data)
}

# Step 1: Apply TPM normalization
count_norm <- normalize_tpm(count_matrix)

# Step 2: Average replicates for each batch
data_avg <- count_norm %>%
  mutate(
    WT_G1 = rowMeans(select(., `WT-1-G1`, `WT-2-G1`)),
    WT_ES = rowMeans(select(., `WT-1-ES`, `WT-2-ES`)),
    WT_MS = `WT-1-MS`,
    WT_LS = `WT-1-LS`,
    
    sol1_5_G1 = rowMeans(select(., `sol1_5-1-G1`, `sol1_5-2-G1`)),
    sol1_5_ES = rowMeans(select(., `sol1_5-1-ES`, `sol1_5-2-ES`)),
    sol1_5_MS = `sol1_5-1-MS`,
    sol1_5_LS = `sol1_5-1-LS`,
    
    sol1_8_G1 = rowMeans(select(., `sol1_8-1-G1`, `sol1_8-2-G1`)),
    sol1_8_ES = rowMeans(select(., `sol1_8-1-ES`, `sol1_8-2-ES`)),
    sol1_8_MS = `sol1_8-1-MS`,
    sol1_8_LS = `sol1_8-1-LS`,
    
    tcx2_1_G1 = rowMeans(select(., `tcx2_1-1-G1`, `tcx2_1-2-G1`)),
    tcx2_1_ES = rowMeans(select(., `tcx2_1-1-ES`, `tcx2_1-2-ES`)),
    tcx2_1_MS = `tcx2_1-1-MS`,
    tcx2_1_LS = `tcx2_1-1-LS`,
    
    tcx2_3_G1 = `tcx2_3-1-G1`,
    tcx2_3_ES = `tcx2_3-1-ES`,
    tcx2_3_MS = `tcx2_3-1-MS`,
    tcx2_3_LS = `tcx2_3-1-LS`
  )

# Step 3: Normalize S-phase by G1 per batch
batches <- c("WT", "sol1_5", "sol1_8", "tcx2_1", "tcx2_3")
for (batch in batches) {
  g1 <- paste0(batch, "_G1")
  for (phase in c("ES", "MS", "LS")) {
    raw <- paste0(batch, "_", phase)
    norm <- paste0(raw, "_norm")
    data_avg[[norm]] <- ifelse(data_avg[[g1]] == 0, 0,
                               ifelse(data_avg[[raw]] / data_avg[[g1]] >= 1,
                                      data_avg[[raw]] / data_avg[[g1]], 0))
  }
}

# Step 4: Classify per batch
assign_phase <- function(es, ms, ls) {
  signals <- c(ES = es, MS = ms, LS = ls)
  if (all(signals == 0)) {
    "Non-replication"
  } else {
    top <- names(signals)[signals >= 0.9 * max(signals)]
    paste(top, collapse = "")
  }
}

for (batch in batches) {
  es <- data_avg[[paste0(batch, "_ES_norm")]]
  ms <- data_avg[[paste0(batch, "_MS_norm")]]
  ls <- data_avg[[paste0(batch, "_LS_norm")]]
  data_avg[[paste0(batch, "_phase")]] <- mapply(assign_phase, es, ms, ls)
}

# Step 5: Final classification by voting
phase_cols <- paste0(batches, "_phase")
data_avg <- data_avg %>%
  rowwise() %>%
  mutate(
    final_phase = {
      phases <- c_across(all_of(phase_cols))
      if (all(phases == "Non-replication")) {
        "Non-replication"
      } else {
        known <- phases[phases != "Non-replication"]
        if (length(unique(known)) == length(known)) {
          "unknown"
        } else {
          names(sort(table(known), decreasing = TRUE))[1]
        }
      }
    }
  ) %>%
  ungroup()

# Step 6: Export GFF3 and CSV
classified <- data_avg %>%
  filter(!final_phase %in% c("Non-replication", "unknown")) %>%
  mutate(
    start = format(start, scientific = FALSE, trim = TRUE),
    end = format(end, scientific = FALSE, trim = TRUE),
    final_phase = gsub("S", "", final_phase),
    annotation = case_when(
      final_phase == "E" ~ "Name=ES;color=#2250F1;",
      final_phase == "EM" ~ "Name=ESMS;color=#28C5CC;",
      final_phase == "M" ~ "Name=MS;color=#1A8A12;",
      final_phase == "ML" ~ "Name=MSLS;color=#FFFD33;",
      final_phase == "L" ~ "Name=LS;color=#FB0018;",
      final_phase == "EL" ~ "Name=ESLS;color=#EA3CF2;",
      final_phase == "EML" ~ "Name=ESMSLS;color=#FAB427;",
      TRUE ~ "Name=unknown;color=#AAAAAA;"
    )
  )

gtf <- classified %>%
  mutate(source = ".", feature = "peaks", score = ".", strand = ".", frame = ".") %>%
  select(chr, source, feature, start, end, score, strand, frame, annotation)

write.table(gtf, "repli_phases.gff3", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.csv(classified, "repli_phase_classified.csv", row.names = FALSE)