############################################################
## 0️⃣ 加载包
############################################################
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(dplyr)
library(stringr)

############################################################
## 1️⃣ 读入 count matrix
############################################################
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

############################################################
## 2️⃣ 计算 TPM
############################################################
length_kb <- gene_length / 1000
rpk <- sweep(count_matrix, 1, length_kb, "/")
tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6

############################################################
## 3️⃣ 定义 sample → genotype
############################################################
sample_group <- data.frame(
  sample = colnames(tpm),
  genotype = c(
    "TCX2", "WT", "WT",
    "SOL1", "SOL1",
    "SOL1", "SOL1",
    "TCX2", "TCX2",
    "TCX2"
  ),
  stringsAsFactors = FALSE
)

############################################################
## 4️⃣ 在重复之间取平均 TPM
############################################################
tpm_mean <- sapply(
  split(sample_group$sample, sample_group$genotype),
  function(samples) {
    rowMeans(tpm[, samples, drop = FALSE])
  }
)
# 行 = gene (MSU)，列 = WT / TCX2 / SOL1

############################################################
## 5️⃣ 计算 mutant / WT 的 log2 fold change
############################################################
eps <- 1
log2fc_mat <- data.frame(
  TCX2_vs_WT = log2((tpm_mean[, "TCX2"] + eps) / (tpm_mean[, "WT"] + eps)),
  SOL1_vs_WT = log2((tpm_mean[, "SOL1"] + eps) / (tpm_mean[, "WT"] + eps))
)

############################################################
## 6️⃣ 读取表观遗传基因注释表（MSU → symbol → type）
############################################################
gene_anno <- read_excel("histone_gene.xlsx") %>%
  dplyr::select(msu, symbol, type) %>%
  dplyr::rename(
    MSU    = msu,
    Symbol = symbol,
    Type   = type
  ) %>%
  filter(!is.na(MSU), !is.na(Symbol)) %>%
  mutate(
    Type   = str_trim(Type),
    Symbol = str_trim(Symbol)
  ) %>%
  distinct()

############################################################
## 7️⃣ 构建 heatmap 矩阵（先用 MSU 对齐）
############################################################
heatmap_mat <- log2fc_mat[
  rownames(log2fc_mat) %in% gene_anno$MSU,
]

# 行 = mutant，列 = gene（MSU）
heatmap_mat <- t(heatmap_mat)

############################################################
## 8️⃣ 将 heatmap 列名从 MSU 映射为 Symbol
############################################################
msu_to_symbol <- gene_anno$Symbol
names(msu_to_symbol) <- gene_anno$MSU

colnames(heatmap_mat) <- make.unique(
  msu_to_symbol[colnames(heatmap_mat)]
)

############################################################
## 9️⃣ 构建 annotation_col（⚠️ rownames = Symbol）
############################################################
annotation_col <- gene_anno %>%
  filter(Symbol %in% colnames(heatmap_mat)) %>%
  mutate(
    Type = factor(
      Type,
      levels = c(
        "DNA Methylation",
        "Histone  Demethylation",
        "Histone  Methylation",
        "Histone Acetyltransferase",
        "siRNA"
      ),
      ordered = TRUE
    )
  ) %>%
  arrange(Type) %>%
  dplyr::select(Symbol, Type) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Symbol
annotation_col$Symbol <- NULL   # ⚠️ 只剩一列：Type

############################################################
## 🔟 按 annotation 顺序同步重排 heatmap
############################################################
heatmap_mat <- heatmap_mat[, rownames(annotation_col), drop = FALSE]

############################################################
## 1️⃣1️⃣ 定义 annotation 颜色（名字必须 = Type）
############################################################
annotation_colors <- list(
  Type = c(
    "DNA Methylation"            = "#4F81BD",
    "Histone  Demethylation"     = "#C0504D",
    "Histone  Methylation"       = "#9BBB59",
    "Histone Acetyltransferase"  = "#8064A2",
    "siRNA"                      = "#F79646"
  )
)

############################################################
## 1️⃣2️⃣ 画 heatmap（Type 现在一定会显示）
############################################################
pdf("Epigenetic_regulators_log2FC_mutant_vs_WT_symbol_withType.pdf",
    width = 12, height = 3)

color_palette <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(11, "RdBu"))
)(100)

pheatmap(
  heatmap_mat,
  color = color_palette,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,   # TCX2_vs_WT / SOL1_vs_WT
  show_colnames = TRUE,   # Symbol
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "log2 Fold Change (mutant / WT)",
  cellwidth  = 10,
  cellheight = 10,
  fontsize_row = 10,
  fontsize_col = 7,
  scale = "none",
  breaks = seq(-2, 2, length.out = 101)
)

dev.off()