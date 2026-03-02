# ===============================
# 1️⃣ 加载包
# ===============================
library(pheatmap)
library(RColorBrewer)

# ===============================
# 2️⃣ 转置 heatmap（行 = mutant，列 = gene）
# ===============================
heatmap_mat_t <- t(heatmap_mat)

# ===============================
# 3️⃣ 构建 gene → cell cycle phase 映射
# ===============================
genes_raw <- c(
  "LOC_Os07g47180","LOC_Os03g50860","LOC_Os10g42950",
  "LOC_Os01g10690","LOC_Os01g36390","LOC_Os05g19270",
  "LOC_Os08g42600","LOC_Os01g67160","LOC_Os01g67160",
  "LOC_Os03g49750","LOC_Os01g13260","LOC_Os01g59120",
  "LOC_Os01g67160","LOC_Os05g41390","LOC_Os06g51110",
  "LOC_Os06g51110","LOC_Os01g63710","LOC_Os01g63710",
  "LOC_Os02g04240","LOC_Os03g41100","LOC_Os05g41080",
  "LOC_Os10g33310"
)

phases_raw <- c(
  "G0","G1","G1","G1/S","G1/S","G1/S","G1/S",
  "G2","G2","G2",
  "M","M","M","M","M","M",
  "S","S","S","S","S","S"
)

gene_phase_df <- data.frame(
  Gene  = genes_raw,
  Phase = phases_raw,
  stringsAsFactors = FALSE
)

# 去重，保留第一次出现
gene_phase_df <- gene_phase_df[!duplicated(gene_phase_df$Gene), ]
rownames(gene_phase_df) <- gene_phase_df$Gene
gene_phase_df$Gene <- NULL

# ===============================
# 4️⃣ 定义 cell cycle 的“生物学顺序”
# ===============================
phase_levels <- c("G0", "G1", "G1/S", "S", "G2", "M")

gene_phase_df$Phase <- factor(
  gene_phase_df$Phase,
  levels = phase_levels,
  ordered = TRUE
)

# ===============================
# 5️⃣ 对齐 heatmap 中的基因，并按 phase 排序
# ===============================
annotation_col <- gene_phase_df[colnames(heatmap_mat_t), , drop = FALSE]

# 按 Phase 排序列（gene）
gene_order <- order(annotation_col$Phase)
heatmap_mat_t <- heatmap_mat_t[, gene_order, drop = FALSE]
annotation_col <- annotation_col[gene_order, , drop = FALSE]

# ===============================
# 6️⃣ 定义 phase 对应的颜色（与你的顺序一致）
# ===============================
annotation_colors <- list(
  Phase = c(
    "G0"   = "#BFBFBF",
    "G1"   = "#4F81BD",
    "G1/S" = "#5FA777",
    "S"    = "#8E83B7",
    "G2"   = "#E09B3D",
    "M"    = "#C84C4C"
  )
)

# ===============================
# 7️⃣ 画 heatmap（已排序）
# ===============================
pdf("heatmap_log2FC_mutants_vs_WT_cellcycle_sorted.pdf",
    width = 12, height = 3)

color_palette <- colorRampPalette(colors = RColorBrewer::brewer.pal(11,"RdBu"))(100)
color_palette <- color_palette %>% rev()

pheatmap(
  heatmap_mat_t,
  color = color_palette,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  show_rownames = TRUE,   # mutant
  show_colnames = TRUE,   # gene
  
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  
  main = "log2 Fold Change vs WT (Cell-cycle ordered)",
  
  cellwidth  = 20,
  cellheight = 20,
  
  fontsize_row = 10,
  fontsize_col = 8,
  
  scale = "none",
  breaks = seq(-2, 2, length.out = 101)
)

dev.off()