library(edgeR)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
set.seed(2025)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/")

# 1) 基因集
cellcycle_genes <- read.csv("DNA_repair_genes.csv")$msu
cellcycle_genes <- unique(trimws(cellcycle_genes))

# 2) 读取和处理
cts <- read.table("./gene_counts.txt", header = TRUE, check.names = FALSE)
mat <- cts[, 7:ncol(cts)]
rownames(mat) <- trimws(cts$Geneid)

group_all <- c("TCX2-3-KO","WT","WT","SOL1-5-KO","SOL1-5-KO",
               "SOL1-8-KO","SOL1-8-KO","TCX2-1-KO","TCX2-1-KO","TCX2-3-KO")

# 3) edgeR
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

# 4) 排序基因列表（降序）
gene_list <- deg_full$logFC
names(gene_list) <- rownames(deg_full)
gene_list <- sort(gene_list, decreasing = TRUE) 

# 5) GSEA
gsea_res <- GSEA(gene_list,
                 TERM2GENE = data.frame(
                   term = "Cell_Cycle",
                   gene = cellcycle_genes
                 ),
                 eps = 0,
                 minGSSize = 2,
                 pvalueCutoff = 1,
                 nPermSimple = 10000)

# 6) 绘图
pdf("GSEA_DNA_Repair.pdf", width = 6, height = 6)
gseaplot2(
  gsea_res,
  geneSetID = 1,
  subplots = 1:3,
  rel_heights = c(1.5, 0.5, 1.0),
  color = "#E64B35FF",
  ES_geom = "line",
  title = "DNA Repair Gene Set Enrichment"
)
dev.off()
