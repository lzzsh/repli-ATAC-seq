setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/")
anno <- read_xlsx("../../ATAC-seq-CR/diffbind/gene_function.xlsx")
anno <- anno[,c(1,10)]
DAR_sol1 <- read.csv("./sol-8_CR_WT_ES_edger_sig_with_cuttag_nearTSS3kb_sorted.csv")
DAR_tcx2 <- read.csv("./TCX2-3_CR_WT_ES_edger_sig_with_cuttag_nearTSS3kb_sorted.csv")
DEG_sol1 <- read.csv("../../RNA-seq/star_all_rawdata/featureCounts_out/edgeR_SOL1-8-KO_vs_WT_sig.csv")
DEG_tcx2 <- read.csv("../../RNA-seq/star_all_rawdata/featureCounts_out/edgeR_TCX2-3-KO_vs_WT_sig.csv")

DAR_sol1_anno <- DAR_sol1 %>%
  left_join(anno, by = c("geneId" = "MSU"))

DAR_tcx2_anno <- DAR_tcx2 %>%
  left_join(anno, by = c("geneId" = "MSU"))

DEG_sol1_anno <- DEG_sol1 %>%
  left_join(anno, by = c("X" = "MSU"))

DEG_tcx2_anno <- DEG_tcx2 %>%
  left_join(anno, by = c("X" = "MSU"))

# 输出 DAR_sol1 注释表
write.csv(DAR_sol1_anno,
          file = "DAR_sol1_anno.csv",
          row.names = FALSE)

# 输出 DAR_tcx2 注释表
write.csv(DAR_tcx2_anno,
          file = "DAR_tcx2_anno.csv",
          row.names = FALSE)

# 输出 DEG_sol1 注释表
write.csv(DEG_sol1_anno,
          file = "DEG_sol1_anno.csv",
          row.names = FALSE)

# 输出 DEG_tcx2 注释表
write.csv(DEG_tcx2_anno,
          file = "DEG_tcx2_anno.csv",
          row.names = FALSE)
