library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(gridExtra)
library(limma)

# 设置工作目录
setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation/")

# 设置复制时期（自动控制 motif 筛选与输出命名）
stage <- "LS"  # 可选 ES / MS / LS

# 读取原始数据
count_matrix <- read.table("./repli_peaks_quan_peak.txt")
colnames(count_matrix) <- c(
  "chr","start","end",
  "WT-1-G1", "WT-2-G1", "WT-1-ES", "WT-2-ES", "WT-1-MS", "WT-1-LS",
  "sol1_5-1-G1", "sol1_5-2-G1", "sol1_5-1-ES", "sol1_5-2-ES", "sol1_5-1-MS", "sol1_5-1-LS",
  "sol1_8-1-G1", "sol1_8-2-G1", "sol1_8-1-ES", "sol1_8-2-ES", "sol1_8-1-MS", "sol1_8-1-LS",
  "tcx2_1-1-G1", "tcx2_1-2-G1", "tcx2_1-1-ES", "tcx2_1-2-ES", "tcx2_1-1-MS", "tcx2_1-1-LS",
  "tcx2_3-1-G1", "tcx2_3-1-ES", "tcx2_3-1-MS", "tcx2_3-1-LS"
)

# TPM 标准化函数
normalize_tpm <- function(data) {
  data$RegionLength <- data$end - data$start + 1
  data[, 4:31] <- lapply(data[, 4:31], function(col) {
    rpk <- col / (data$RegionLength / 1000)
    scaling_factor <- sum(rpk) / 1e6
    tpm <- rpk / scaling_factor
    return(tpm)
  })
  return(data)
}
count_matrix <- normalize_tpm(count_matrix)

# G1 normalization（按组）
all_cols <- colnames(count_matrix)[4:31]
groups <- unique(sub("-[A-Z]+[12]?$", "", all_cols))  # eg: WT, sol1_5
for (grp in groups) {
  g1_cols <- grep(paste0("^", grp, "-G1"), colnames(count_matrix), value = TRUE)
  g1_mean <- rowMeans(count_matrix[, g1_cols, drop = FALSE])
  for (s in c("ES", "MS", "LS")) {
    stage_cols <- grep(paste0("^", grp, "-", s), colnames(count_matrix), value = TRUE)
    count_matrix[, stage_cols] <- count_matrix[, stage_cols] / (g1_mean + 1e-6)
  }
}

rownames(count_matrix) <- paste(count_matrix$chr, count_matrix$start, count_matrix$end)

# 读入 motif 并筛选指定复制阶段
sol1 <- read.table("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/fimo/meme/SOL1_overlap_peak.bed") %>%
  filter(V4 == stage)
colnames(sol1)[1:3] <- c("chr","start","end")

tcx2 <- read.table("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/fimo/meme/TCX2_overlap_peak.bed") %>%
  filter(V4 == stage)
colnames(tcx2)[1:3] <- c("chr","start","end")

# 添加 motif 标记列
count_matrix <- count_matrix %>%
  mutate(
    motif_SOL1 = as.integer(paste(chr, start, end) %in% paste(sol1$chr, sol1$start, sol1$end)),
    motif_TCX2 = as.integer(paste(chr, start, end) %in% paste(tcx2$chr, tcx2$start, tcx2$end))
  )

# log2 + quantile normalization
expr_matrix <- as.matrix(count_matrix[, -(1:3)])
log2_matrix <- log2(expr_matrix + 1)
qn_matrix <- normalizeBetweenArrays(log2_matrix, method = "quantile")
qn_df <- as.data.frame(qn_matrix)

# 添加 motif 信息列
qn_df$motif_SOL1 <- count_matrix$motif_SOL1
qn_df$motif_TCX2 <- count_matrix$motif_TCX2

# 平均表达值（仅 ES 阶段）
qn_df$WT_ES     <- rowMeans(qn_df[, c("WT-1-ES", "WT-2-ES")])
qn_df$sol1_5_ES <- rowMeans(qn_df[, c("sol1_5-1-ES", "sol1_5-2-ES")])
qn_df$sol1_8_ES <- rowMeans(qn_df[, c("sol1_8-1-ES", "sol1_8-2-ES")])
qn_df$tcx2_1_ES <- rowMeans(qn_df[, c("tcx2_1-1-ES", "tcx2_1-2-ES")])
qn_df$tcx2_3_ES <- qn_df[, "tcx2_3-1-ES"]

# 绘图准备
mutant_cols <- c("sol1_5_ES", "sol1_8_ES", "tcx2_1_ES", "tcx2_3_ES")
motif_list <- list(SOL1 = "motif_SOL1", TCX2 = "motif_TCX2", All = NA)

# 绘图函数
plot_scatter_all <- function(motif_name, motif_col = NULL) {
  plots <- list()
  ratios <- c()
  
  for (mutant in mutant_cols) {
    if (!is.na(motif_col)) {
      df_plot <- qn_df %>%
        select(WT = WT_ES, MUT = all_of(mutant), motif_flag = all_of(motif_col)) %>%
        mutate(motif = ifelse(motif_flag == 1, motif_name, "Other")) %>%
        mutate(motif = factor(motif, levels = c("Other", motif_name)))
      
      df_motif <- df_plot %>% filter(motif == motif_name)
    } else {
      df_plot <- qn_df %>%
        select(WT = WT_ES, MUT = all_of(mutant)) %>%
        mutate(motif = "All") %>%
        mutate(motif = factor(motif, levels = "All"))
      
      df_motif <- df_plot
    }
    
    up_ratio <- mean(df_motif$MUT > df_motif$WT, na.rm = TRUE)
    down_ratio <- 1 - up_ratio
    ratios <- c(ratios, up_ratio)
    
    min_val <- min(c(df_plot$WT, df_plot$MUT), na.rm = TRUE)
    max_val <- max(c(df_plot$WT, df_plot$MUT), na.rm = TRUE)
    
    p <- ggplot() +
      geom_point(data = df_plot %>% filter(motif == "Other"), 
                 aes(x = WT, y = MUT), color = "gray80", alpha = 0.4, size = 0.6) +
      geom_point(data = df_plot %>% filter(motif == motif_name), 
                 aes(x = WT, y = MUT), color = "steelblue", alpha = 0.9, size = 1.4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      coord_fixed(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
      labs(
        title = paste0(motif_name, " motif: WT vs ", mutant),
        subtitle = sprintf("motif only: n = %d | Up: %.1f%% | Down: %.1f%%", 
                           nrow(df_motif), up_ratio * 100, down_ratio * 100),
        x = "WT ES (avg)", y = paste0(mutant, " (avg)")
      ) +
      theme_minimal(base_size = 12)
    
    plots[[mutant]] <- p
  }
  
  message(paste0("[", motif_name, "] Proportion of points above y=x (motif only):"))
  print(setNames(round(ratios, 3), mutant_cols))
  
  return(plots)
}

# 输出所有 PDF 图
for (motif in names(motif_list)) {
  motif_col <- motif_list[[motif]]
  p_list <- plot_scatter_all(motif, motif_col)
  pdf(paste0("scatter_all_", motif, "_motif_", stage, "_G1norm.pdf"), width = 11, height = 8)
  grid.arrange(grobs = p_list, ncol = 2)
  dev.off()
}