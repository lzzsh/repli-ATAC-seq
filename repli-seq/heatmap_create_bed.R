library(dplyr)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/heatmap/")

all <- read.table("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT_all_org.gff3")
sol1 <- read.table("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/fimo/meme/SOL1_overlap_peak.bed")
tcx2 <- read.table("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/fimo/meme/TCX2_overlap_peak.bed")

save_stage_regions <- function(df, label, stages = c("ES", "MS", "LS")) {
  for (stage in stages) {
    subset_df <- df[df$V4 == stage, 1:3]
    filename <- paste0(label, "_", stage, ".bed")
    write.table(subset_df, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}

save_stage_regions(sol1, label = "sol1")
save_stage_regions(tcx2, label = "tcx2")
save_stage_regions(all, label = "all")
