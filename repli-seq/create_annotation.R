setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation")

library(ggalluvial)
library(ggplot2)
library(dplyr)
library(tidyr)

df <- read.csv("repli_phase_classified.csv")

write_phase_gff3 <- function(df, sample, output_dir = ".") {
  phase_col <- paste0(sample, "_phase")
  out_df <- df %>%
    filter(!is.na(.data[[phase_col]]), .data[[phase_col]] != "Non-replication", .data[[phase_col]] != "unknown") %>%
    mutate(
      start = ifelse(start == 0, 1, start),  # start=0 改为 1
      phase_clean = gsub("S", "", .data[[phase_col]])
    ) %>%
    mutate(annotation = case_when(
      phase_clean == "E" ~ "Name=ES;color=#2250F1;",
      phase_clean == "EM" ~ "Name=ESMS;color=#28C5CC;",
      phase_clean == "M" ~ "Name=MS;color=#1A8A12;",
      phase_clean == "ML" ~ "Name=MSLS;color=#FFFD33;",
      phase_clean == "L" ~ "Name=LS;color=#FB0018;",
      phase_clean == "EL" ~ "Name=ESLS;color=#EA3CF2;",
      phase_clean == "EML" ~ "Name=ESMSLS;color=#FAB427;",
      TRUE ~ "Name=Unknown;color=#000000;"
    )) %>%
    mutate(
      source = ".",
      feature = "peaks",
      score = ".",
      strand = ".",
      frame = "."
    ) %>%
    select(chr, source, feature, start, end, score, strand, frame, annotation)
  
  output_file <- file.path(output_dir, paste0(sample, "_replication_phase.gff3"))
  write.table(out_df, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message("Saved: ", output_file)
}

# 样本列表
samples <- c("WT", "sol1_5", "sol1_8", "tcx2_1", "tcx2_3")

# 生成每个样本的 GFF3 文件
for (s in samples) {
  write_phase_gff3(df, s)
}