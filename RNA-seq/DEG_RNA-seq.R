############################################################
## edgeR DEG + Heatmap (保持原版计算不变 + ComplexHeatmap anno_mark 引导线)
## 说明：
## 1) edgeR计算、DEG筛选、DEG导出：完全沿用你原版（不改）
## 2) Heatmap用 FDR<0.05 的基因 + 强制包含 my_genes
## 3) 右侧用 anno_mark 给 my_genes 画引导线标注
## 4) 标注顺序严格按 show_genes（my_genes）给定顺序
############################################################

suppressPackageStartupMessages({
  library(edgeR)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/")

#------------------------------------------------------------
# 1. 读入 count matrix（完全沿用你原版）
#------------------------------------------------------------
counts <- read.table("gene_counts.txt", header = TRUE)

count_matrix <- counts[, 7:ncol(counts)]
rownames(count_matrix) <- counts$Geneid
gene_length <- counts$Length

colnames(count_matrix) <- c(
  "TCX2-3-KO.1", "WT.1", "WT.2",
  "SOL1-5-KO.1", "SOL1-5-KO.2",
  "SOL1-8-KO.1", "SOL1-8-KO.2",
  "TCX2-1-KO.1", "TCX2-1-KO.2",
  "TCX2-3-KO.2"
)

# TPM（可选，用于后续展示）
length_kb <- gene_length / 1000
rpk <- sweep(count_matrix, 1, length_kb, "/")
tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6
write.csv(tpm, "./gene_tpm_matrix.csv", quote = FALSE)

group_all <- c(
  "TCX2-3-KO", "WT", "WT",
  "SOL1-5-KO", "SOL1-5-KO",
  "SOL1-8-KO", "SOL1-8-KO",
  "TCX2-1-KO", "TCX2-1-KO",
  "TCX2-3-KO"
)

#------------------------------------------------------------
# 2. edgeR + Heatmap 主函数（edgeR部分完全照抄你原版）
#------------------------------------------------------------
run_edgeR <- function(count_matrix, group_all, target, ref = "WT", show_genes = NULL) {
  
  # 选择样本（原版）
  sel <- group_all %in% c(ref, target)
  count_sub <- count_matrix[, sel]
  group <- factor(group_all[sel], levels = c(ref, target))
  
  # edgeR 标准流程（原版）
  y <- DGEList(counts = count_sub, group = group)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ group)
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  
  tab <- topTags(lrt, n = Inf)$table
  
  #----------------------------------------------------------
  # (1) DEG 导出（严格阈值）（原版）
  #----------------------------------------------------------
  deg <- tab %>% filter(FDR < 0.05 & abs(logFC) > 1)
  
  ## 这里保持你原版 write.csv 行为：写出时保留行名=GeneID
  write.csv(deg,
            paste0("edgeR_", target, "_vs_", ref, "_sig.csv"),
            quote = FALSE)
  
  #----------------------------------------------------------
  # (2) Heatmap 行基因：FDR < 0.05 + 强制包含 show_genes（原版）
  #----------------------------------------------------------
  deg_fdr <- tab %>% filter(FDR < 0.05)
  genes_for_heatmap <- rownames(deg_fdr)
  
  if (!is.null(show_genes)) {
    genes_for_heatmap <- union(genes_for_heatmap, show_genes)
  }
  
  logCPM <- cpm(y, log = TRUE, prior.count = 1)
  sig_expr <- logCPM[rownames(logCPM) %in% genes_for_heatmap, , drop = FALSE]
  
  if (nrow(sig_expr) < 2) {
    message(target, ": not enough genes for heatmap")
    return(nrow(deg))
  }
  
  # Z-score（原版）
  zscore <- t(scale(t(sig_expr)))
  zscore <- zscore[apply(zscore, 1, function(x) all(is.finite(x))), , drop = FALSE]
  
  #----------------------------------------------------------
  # (3) 绘图：ComplexHeatmap + anno_mark（只改绘图，不改计算）
  #     关键：标注顺序严格按 show_genes 原始顺序
  #----------------------------------------------------------
  
  # 列注释（等价于你原来的 ann_col + ann_colors）
  ha_col <- HeatmapAnnotation(
    Condition = group,
    col = list(Condition = c(
      WT = "#3498db",
      setNames("#e74c3c", target)
    )),
    annotation_name_side = "left"
  )
  
  # 颜色（尽量贴近你原来 pheatmap 的 6 色渐变）
  color_palette <- colorRampPalette(
    rev(c("#D93F49", "#E28187", "#EBBFC2",
          "#D5E1E3", "#AFC9CF", "#8FB4BE"))
  )(100)
  
  # 用 zscore 的范围生成 colorRamp2（更接近 pheatmap 的效果）
  brks <- seq(min(zscore, na.rm = TRUE), max(zscore, na.rm = TRUE), length.out = 100)
  col_fun <- colorRamp2(brks, color_palette)
  
  # ---- 计算需要标注的行位置（按 show_genes 顺序）----
  row_anno <- NULL
  if (!is.null(show_genes)) {
    mark_genes <- show_genes[show_genes %in% rownames(zscore)]  # 保持原顺序
    if (length(mark_genes) > 0) {
      mark_at <- match(mark_genes, rownames(zscore))           # 对应行号（按 show_genes 顺序）
      row_anno <- rowAnnotation(
        mark = anno_mark(
          at = mark_at,
          labels = mark_genes,                 # ✅ 顺序=show_genes
          labels_gp = gpar(fontsize = 10),
          link_gp = gpar(lwd = 1),
          padding = unit(1, "mm")
        )
      )
    }
  }
  
  # 输出 PDF
  pdf(paste0("heatmap_", target, "_vs_", ref, "_markline.pdf"),
      width = 8, height = 10)
  
  ht <- Heatmap(
    zscore,
    name = "Z",
    col = col_fun,
    top_annotation = ha_col,
    show_row_names = FALSE,          # 关闭默认行名，避免拥挤
    show_column_names = TRUE,
    column_names_rot = 90,
    
    # 让聚类方式更接近你原来的 pheatmap 参数
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    
    row_title = paste0(target, " vs ", ref),
    column_title = paste0(target, " vs ", ref, "\nDifferentially Expressed Genes (FDR < 0.05)")
  )
  
  if (!is.null(row_anno)) {
    draw(ht + row_anno, heatmap_legend_side = "right", annotation_legend_side = "right")
  } else {
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  }
  
  dev.off()
  
  return(nrow(deg))  # ✅ DEG数量=你原版
}

#------------------------------------------------------------
# 3. 运行所有比较
#------------------------------------------------------------
my_genes <- c(
  "LOC_Os12g31370",
  "LOC_Os12g04980",
  "LOC_Os05g41880",
  "LOC_Os04g41110",
  "LOC_Os05g50410",
  "LOC_Os02g56130"
)

n1 <- run_edgeR(count_matrix, group_all, target = "SOL1-5-KO", show_genes = my_genes)
n2 <- run_edgeR(count_matrix, group_all, target = "SOL1-8-KO", show_genes = my_genes)
n3 <- run_edgeR(count_matrix, group_all, target = "TCX2-1-KO", show_genes = my_genes)
n4 <- run_edgeR(count_matrix, group_all, target = "TCX2-3-KO", show_genes = my_genes)

summary <- data.frame(
  Comparison = c("SOL1-5-KO_vs_WT", "SOL1-8-KO_vs_WT", "TCX2-1-KO_vs_WT", "TCX2-3-KO_vs_WT"),
  DEG_count = c(n1, n2, n3, n4)
)
write.csv(summary, "edgeR_DEG_summary.csv", row.names = FALSE, quote = FALSE)