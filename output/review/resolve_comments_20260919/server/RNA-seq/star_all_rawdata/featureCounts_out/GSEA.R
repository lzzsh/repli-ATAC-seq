library(edgeR)
library(fgsea)
library(tidyverse)
set.seed(2025)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/")

cellcycle_genes <- read.csv("Strictly_Matched_Cell_Cycle_Genes.csv")$msu
gene_sets <- list(Cell_Cycle = unique(trimws(cellcycle_genes)))

cts <- read.table("./gene_counts.txt", header = TRUE, check.names = FALSE)
mat <- cts[, 7:ncol(cts)]
rownames(mat) <- trimws(cts$Geneid)

group_all <- c("TCX2-3-KO","WT","WT","SOL1-5-KO","SOL1-5-KO",
               "SOL1-8-KO","SOL1-8-KO","TCX2-1-KO","TCX2-1-KO","TCX2-3-KO")

sel <- group_all %in% c("WT","TCX2-3-KO")
y  <- DGEList(counts = mat[, sel], group = factor(group_all[sel], levels=c("WT","TCX2-3-KO")))
keep <- filterByExpr(y); y <- y[keep,, keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ y$samples$group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
deg_full <- topTags(lrt, n = Inf)$table
deg_full <- deg_full[!duplicated(rownames(deg_full)), ]
rownames(deg_full) <- trimws(rownames(deg_full))

gene_list <- deg_full$logFC
names(gene_list) <- rownames(deg_full)
gene_list <- sort(gene_list, decreasing = FALSE)

intersect_genes <- intersect(names(gene_list), gene_sets$Cell_Cycle)
gsea_res <- fgsea(pathways = gene_sets, stats = gene_list, nperm = 10000)

# ===== 改进的可视化 =====
enrich_plot <- plotEnrichment(gene_sets$Cell_Cycle, gene_list)
pd <- ggplot_build(enrich_plot)$data[[1]]

# 美化版本
ggplot(pd, aes(x = x, y = y)) +
  # 背景
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "#FAFAFA", color = NA) +
  
  # 富集曲线
  geom_line(color = "#E64B35FF", linewidth = 1.1) +
  
  # 基因位置标记（下方rug plot）
  geom_rug(aes(x = x), sides = "b", color = "#E64B35FF",
           alpha = 0.6, length = unit(0.02, "npc")) +
  
  # 零线
  geom_hline(yintercept = 0, color = "#333333", linewidth = 0.5, linetype = "dashed") +
  
  # 添加阴影区域（可选，突出富集）
  geom_area(fill = "#E64B35FF", alpha = 0.15) +
  
  # 标签和主题
  labs(
    title = "Cell Cycle Gene Set Enrichment",
    subtitle = paste0("NES = ", round(gsea_res$NES, 2),
                      " | FDR = ", signif(gsea_res$padj, 2),
                      " | Size = ", gsea_res$size),
    x = "Gene Rank (WT → KO)",
    y = "Enrichment Score",
    caption = paste0("Matched genes: ", length(intersect_genes))
  ) +
  
  # 美观主题
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#555555", margin = margin(b = 10)),
    plot.caption = element_text(size = 10, color = "#888888", hjust = 1),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10, color = "#333333"),
    panel.border = element_rect(color = "#CCCCCC", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "#F5F5F5", linewidth = 0.2),
    plot.margin = margin(15, 15, 15, 15)
  )
