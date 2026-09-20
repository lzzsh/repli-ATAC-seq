library(clusterProfiler)
library(enrichplot)
library(tidyverse)
set.seed(2025)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/")

# ================================
# 1. 读取 Arabidopsis→Rice 转换表
# ================================
trans <- read.csv("./rice2ara.csv")  # 包含 rice (MSU) 和 tair 两列
trans$tair <- trimws(trans$tair)
trans$rice <- trimws(trans$rice)

# ================================
# 2. 读取 Arabidopsis 的 SOG1 target
# ================================
cellcycle_genes <- read.csv("SOG1_taget.csv", header = FALSE)$V1
cellcycle_genes <- trimws(cellcycle_genes)

# 转换为水稻 ID
cellcycle_rice <- trans %>%
  filter(tair %in% cellcycle_genes) %>%
  pull(rice) %>%
  unique()

cat("找到的水稻基因数量：", length(cellcycle_rice), "\n")

# ================================
# 3. 读取水稻计数矩阵
# ================================
cts <- read.table("./gene_counts.txt", header = TRUE, check.names = FALSE)
mat <- cts[, 7:ncol(cts)]
rownames(mat) <- trimws(cts$Geneid)

group_all <- c("TCX2-3-KO","WT","WT","SOL1-5-KO","SOL1-5-KO",
               "SOL1-8-KO","SOL1-8-KO","TCX2-1-KO","TCX2-1-KO","TCX2-3-KO")

# ================================
# 4. edgeR 差异分析
# ================================
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

# ================================
# 5. 排序基因列表（GSEA 输入）
# ================================
gene_list <- deg_full$logFC
names(gene_list) <- trimws(rownames(deg_full))
gene_list <- sort(gene_list, decreasing = TRUE)

# ================================
# 6. 过滤未出现在 RNA-seq 中的水稻基因
# ================================
cellcycle_rice <- cellcycle_rice[cellcycle_rice %in% names(gene_list)]

cat("用于 GSEA 的水稻基因数量：", length(cellcycle_rice), "\n")

# ================================
# 7. GSEA 分析
# ================================
gsea_res <- GSEA(gene_list,
                 TERM2GENE = data.frame(
                   term = "SOG1_Direct_Targets",
                   gene = cellcycle_rice
                 ),
                 eps = 0,
                 minGSSize = 2,
                 pvalueCutoff = 1,
                 nPermSimple = 10000)

# ================================
# 8. 绘图
# ================================
pdf("GSEA_SOG1_Targets_TCX2_KO.pdf", width = 6, height = 6)
gseaplot2(
  gsea_res,
  geneSetID = 1,
  subplots = 1:3,
  rel_heights = c(1.5, 0.5, 1.0),
  color = "#E64B35FF",
  ES_geom = "line",
  title = "GSEA of SOG1 Direct Targets in TCX2-KO"
)
dev.off()