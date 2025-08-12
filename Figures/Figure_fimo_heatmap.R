## ================= 依赖 =================
req <- c("dplyr","tidyr","ComplexHeatmap","circlize","readr","stats","grDevices","grid")
for (p in req) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(readr)
library(stats)
library(grDevices)
library(grid)
# library(ggplot2) # 删除这行

## ================= 输入 =================
## 1) tf_table：行=TF；列=ES/ESMS/MS/LS；值=事件数（同一peak内不同变体可多次计）
##    若已在环境中就注释掉下面读入；否则按需指定你的文件路径：
tf_table <- read.table("~/Desktop/Rfiles/peak_unit/tf_location_classfied_org.bed",
                       header = TRUE, sep = "\t", row.names = 1)

## 2) 各复制时期总暴露量（推荐用 bp 总长度；无则用区域总数）
##    用你的真实数值替换：
total_bp <- c(
  ES   = 10008822,
  ESMS = 1659300,
  MS   = 14546125,
  LS   = 18878938
)

## 3) TF 家族映射：两列 [TF, TF_family]
##    你给的文件：/Users/lzz/Desktop/Rfiles/fimo.bed（含 motif_alt_id, TF_family）
##    这里把 motif_alt_id 当作 TF 名；若你的 TF 名另有来源，请相应修改。
tf_family <- read.table("/Users/lzz/Desktop/Rfiles/fimo.bed", header = TRUE, sep = "\t", check.names = FALSE)
tf_family <- tf_family[, c(2, 3)]                           # 取 motif_alt_id, TF_family 两列
tf_family <- dplyr::distinct(tf_family, motif_alt_id, TF_family, .keep_all = TRUE)
colnames(tf_family)[1] <- "TF"

## ================= 对齐与检查 =================
stopifnot(is.data.frame(tf_table), nrow(tf_table) > 0)
periods <- colnames(tf_table)
if (!all(periods %in% names(total_bp))) {
  stop("total_bp 未覆盖 tf_table 的所有列，请补齐：缺失 -> ",
       paste(setdiff(periods, names(total_bp)), collapse = ", "))
}
total_bp <- total_bp[periods]         # 对齐顺序
TFs <- rownames(tf_table)

## ================= Poisson 比率检验 =================
pval_mat   <- matrix(NA_real_, nrow = nrow(tf_table), ncol = ncol(tf_table),
                     dimnames = list(TFs, periods))
log2rr_mat <- matrix(NA_real_, nrow = nrow(tf_table), ncol = ncol(tf_table),
                     dimnames = list(TFs, periods))

for (i in seq_len(nrow(tf_table))) {
  total_hits_tf <- sum(tf_table[i, ], na.rm = TRUE)
  for (j in seq_len(ncol(tf_table))) {
    fg_hits <- tf_table[i, j]
    bg_hits <- total_hits_tf - fg_hits
    
    fg_T <- total_bp[j]
    bg_T <- sum(total_bp) - fg_T
    
    if (fg_T <= 0 || bg_T <= 0 || (fg_hits + bg_hits) == 0) {
      pval_mat[i, j]   <- NA_real_
      log2rr_mat[i, j] <- NA_real_
    } else {
      pt <- poisson.test(c(fg_hits, bg_hits), T = c(fg_T, bg_T), alternative = "two.sided")
      pval_mat[i, j] <- pt$p.value
      
      rate_fg <- fg_hits / fg_T
      rate_bg <- bg_hits / bg_T
      log2rr_mat[i, j] <- log2((rate_fg + 1e-12) / (rate_bg + 1e-12))
    }
  }
}

## ================= 显著性与多重校正 =================
logp_mat <- -log10(pval_mat)

## 按TF分别校正（更合理）
q_mat <- pval_mat
for(i in 1:nrow(pval_mat)) {
  q_mat[i, ] <- p.adjust(pval_mat[i, ], method = "BH")
}
logq_mat <- -log10(q_mat)

## 偏好复制时期（以 -log10(q) 最大为准）
preferred_period <- colnames(logq_mat)[apply(logq_mat, 1, function(x) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) < -log10(0.05)) {
    NA_integer_
  } else {
    which.max(x)
  }
})]
pref_tbl <- data.frame(TF = rownames(logq_mat), preferred_period, row.names = NULL)

## ================= 合并 TF 家族并准备分片信息 =================
## 按 logq_mat 的行顺序，给每个 TF 找家族；缺失/空置为 Unknown
fam_vec <- tf_family$TF_family[match(rownames(logq_mat), tf_family$TF)]
fam_vec[is.na(fam_vec) | fam_vec == ""] <- "Unknown"

## （可选）基于家族“平均模式”决定家族切片顺序，使得家族块在热图中更有序
fam_levels <- {
  fam_centroid <- sapply(split.data.frame(as.data.frame(logq_mat), fam_vec),
                         function(m) colMeans(m, na.rm = TRUE))
  hc <- hclust(dist(t(fam_centroid)))
  colnames(fam_centroid)[hc$order]
}
fam_factor <- factor(fam_vec, levels = fam_levels)

## 行侧注释颜色表（家族彩带）
fam_unique <- levels(fam_factor)

## 使用指定的高对比度颜色方案
specified_colors <- c(
  '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
  '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', 
  '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', 
  '#000075', '#a9a9a9', '#ffffff', '#000000', '#ff6b35', '#004e89',
  '#1a936f', '#88d4ab', '#c6426e', '#642ca9', '#f18f01'
)

n_families <- length(fam_unique)

# 如果家族数量超过颜色数量，循环使用
family_colors <- setNames(
  specified_colors[((1:n_families - 1) %% length(specified_colors)) + 1], 
  fam_unique
)

# 特殊处理Unknown家族（如果存在），使用灰色
if ("Unknown" %in% fam_unique) {
  family_colors["Unknown"] <- "#a9a9a9"  # 使用你指定颜色中的灰色
}

# 输出颜色分配信息
cat("TF家族颜色分配:\n")
for(i in 1:length(family_colors)) {
  cat(sprintf("%s: %s\n", names(family_colors)[i], family_colors[i]))
}

row_anno <- rowAnnotation(
  Family = fam_factor,
  col = list(Family = family_colors),
  annotation_legend_param = list(title = "TF family")
)

## ================= 可视化（交换行列） =================
## 固定列顺序（按 ES → ESMS → MS → LS）
desired_order <- c("ES", "ESMS", "MS", "LS")
logq_mat    <- logq_mat[, desired_order]
log2rr_mat  <- log2rr_mat[, desired_order]
q_mat       <- q_mat[, desired_order]

## 转置矩阵：行=period，列=TF
logq_mat_t   <- t(logq_mat)
log2rr_mat_t <- t(log2rr_mat)
q_mat_t      <- t(q_mat)

## 家族分组转移到列方向（列=TF）
col_fam_factor  <- fam_factor
col_family_cols <- family_colors  # 使用相同的颜色方案

## 1) 显著性热图（-log10(q)，旋转版）
col_fun_q <- colorRamp2(
  c(0, stats::quantile(logq_mat_t, 0.5, na.rm = TRUE), max(logq_mat_t, na.rm = TRUE)),
  c("#F7FBFF", "#6BAED6", "#08306B")
)

ht_q <- Heatmap(
  logq_mat_t,
  name = "-log10(q)",
  col = col_fun_q,
  na_col = "grey90",
  cluster_rows = FALSE,                 # 行=period，通常不聚类
  cluster_columns = TRUE,               # 列=TF，可聚类
  column_split = col_fam_factor,        # << 按家族分列分片
  cluster_column_slices = TRUE,
  gap = grid::unit(2, "mm"),
  top_annotation = HeatmapAnnotation(
    Family = col_fam_factor,
    col = list(Family = col_family_cols),
    annotation_legend_param = list(title = "TF family")
  ),
  row_title = "Replication Phase",
  column_title = "TF",
  show_column_names = FALSE
)
draw(ht_q)

## 2) 效应量热图（log2 rate ratio，旋转版）
rng <- max(abs(log2rr_mat_t), na.rm = TRUE)
col_fun_e <- colorRamp2(c(-rng, 0, rng), c("#7F0000", "#FFFFFF", "#00441B"))

ht_e <- Heatmap(
  log2rr_mat_t,
  name = "log2 RR",
  col = col_fun_e,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  column_split = col_fam_factor,
  cluster_column_slices = TRUE,
  gap = grid::unit(2, "mm"),
  top_annotation = HeatmapAnnotation(
    Family = col_fam_factor,
    col = list(Family = col_family_cols),
    annotation_legend_param = list(title = "TF family")
  ),
  row_title = "Replication Phase",
  column_title = "TF",
  show_column_names = FALSE
)
draw(ht_e)

output_dir <- "/Users/lzz/Documents/GitHub/repli-ATAC-seq/output/Figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
pdf(file.path(output_dir, "heatmap_log10Q.pdf"), 
    width = 15, height = 6)
draw(ht_q)
dev.off()

pdf(file.path(output_dir, "heatmap_log2RR.pdf"), 
    width = 15, height = 6)
draw(ht_e)
dev.off()

## ================= 导出结果 =================
write.table(data.frame(TF = rownames(logp_mat), logp_mat),
            file = "TF_phase_log10P.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(TF = rownames(logq_mat), logq_mat),
            file = "TF_phase_log10Q.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(TF = rownames(log2rr_mat), log2rr_mat),
            file = "TF_phase_log2RR.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pref_tbl, file = "TF_preferred_phase.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

message("完成：旋转版热图已生成，并导出 TF_phase_log10P.tsv / TF_phase_log10Q.tsv / TF_phase_log2RR.tsv / TF_preferred_phase.tsv")