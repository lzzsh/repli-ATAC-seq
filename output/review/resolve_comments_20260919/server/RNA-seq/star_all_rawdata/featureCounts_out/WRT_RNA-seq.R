library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/")

# =========================================================
# 1. 读入 RNA TPM matrix + TSS
# =========================================================
tpm_matrix <- fread("./gene_tpm_matrix.csv", check.names = FALSE)
tss_site <- fread(
  "/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/TSS.bed",
  header = FALSE
)

tpm_matrix <- cbind(tss_site[, 1:3], tpm_matrix)
colnames(tpm_matrix)[1:4] <- c("chr", "start", "end", "gene_id")
tpm_matrix[, start := as.numeric(start)]
tpm_matrix[, end   := as.numeric(end)]

# =========================================================
# 2. 读入 repli-seq 1-kb bin 原始信号
# =========================================================
repli_matrix <- fread(
  "/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-seq-CR/segmentation/repli_peaks_quan.txt",
  header = FALSE
)

colnames(repli_matrix) <- c(
  "chr","start","end",
  "WT-1-G1","WT-2-G1","WT-1-ES","WT-2-ES","WT-1-MS","WT-1-LS",
  "sol1_5-1-G1","sol1_5-2-G1","sol1_5-1-ES","sol1_5-2-ES","sol1_5-1-MS","sol1_5-1-LS",
  "sol1_8-1-G1","sol1_8-2-G1","sol1_8-1-ES","sol1_8-2-ES","sol1_8-1-MS","sol1_8-1-LS",
  "tcx2_1-1-G1","tcx2_1-2-G1","tcx2_1-1-ES","tcx2_1-2-ES","tcx2_1-1-MS","tcx2_1-1-LS",
  "tcx2_3-1-G1","tcx2_3-1-ES","tcx2_3-1-MS","tcx2_3-1-LS"
)
repli_matrix[, start := as.numeric(start)]
repli_matrix[, end   := as.numeric(end)]

# =========================================================
# 3. TPM normalization for repli bins
# =========================================================
normalize_tpm <- function(data, value_cols = 4:31) {
  data <- copy(data)
  data[, RegionLength := end - start + 1]
  for (j in value_cols) {
    rpk            <- data[[j]] / (data$RegionLength / 1000)
    scaling_factor <- sum(rpk, na.rm = TRUE) / 1e6
    data[[j]]      <- rpk / scaling_factor
  }
  return(data)
}

count_norm <- normalize_tpm(repli_matrix)

# =========================================================
# 4. 平均 repli 重复（cpp8 = tcx2_3，cpp11 = sol1_8）
# =========================================================
data_avg <- count_norm %>%
  mutate(
    WT_G1     = rowMeans(select(., `WT-1-G1`,     `WT-2-G1`),     na.rm = TRUE),
    WT_ES     = rowMeans(select(., `WT-1-ES`,     `WT-2-ES`),     na.rm = TRUE),
    WT_MS     = `WT-1-MS`,
    WT_LS     = `WT-1-LS`,
    sol1_8_G1 = rowMeans(select(., `sol1_8-1-G1`, `sol1_8-2-G1`), na.rm = TRUE),
    sol1_8_ES = rowMeans(select(., `sol1_8-1-ES`, `sol1_8-2-ES`), na.rm = TRUE),
    sol1_8_MS = `sol1_8-1-MS`,
    sol1_8_LS = `sol1_8-1-LS`,
    tcx2_3_G1 = `tcx2_3-1-G1`,
    tcx2_3_ES = `tcx2_3-1-ES`,
    tcx2_3_MS = `tcx2_3-1-MS`,
    tcx2_3_LS = `tcx2_3-1-LS`
  )

# =========================================================
# 5. S/G1 normalization
# =========================================================
normalize_by_g1 <- function(s, g1, pseudo = 1e-6) s / (g1 + pseudo)

for (batch in c("WT", "sol1_8", "tcx2_3")) {
  data_avg[[paste0(batch, "_ES_norm")]] <- normalize_by_g1(
    data_avg[[paste0(batch, "_ES")]], data_avg[[paste0(batch, "_G1")]])
  data_avg[[paste0(batch, "_MS_norm")]] <- normalize_by_g1(
    data_avg[[paste0(batch, "_MS")]], data_avg[[paste0(batch, "_G1")]])
  data_avg[[paste0(batch, "_LS_norm")]] <- normalize_by_g1(
    data_avg[[paste0(batch, "_LS")]], data_avg[[paste0(batch, "_G1")]])
}

# =========================================================
# 6. 计算 WRT 及 deltaRT
# =========================================================
calc_wrt <- function(es, ms, ls, pseudo = 1e-6) {
  (0.5 * ms + 1 * ls) / (es + ms + ls + pseudo)
}

for (batch in c("WT", "sol1_8", "tcx2_3")) {
  data_avg[[paste0(batch, "_WRT")]] <- calc_wrt(
    data_avg[[paste0(batch, "_ES_norm")]],
    data_avg[[paste0(batch, "_MS_norm")]],
    data_avg[[paste0(batch, "_LS_norm")]]
  )
}

data_avg <- data_avg %>%
  mutate(
    deltaRT_cpp8  = tcx2_3_WRT - WT_WRT,
    deltaRT_cpp11 = sol1_8_WRT - WT_WRT
  )

# =========================================================
# 7. RNA：计算均值与 log2FC（排除 WT_mean == 0）
# =========================================================
expr_dt <- copy(tpm_matrix)

for (cc in c("TCX2-3-KO.1","TCX2-3-KO.2","WT.1","WT.2","SOL1-8-KO.1","SOL1-8-KO.2")) {
  expr_dt[[cc]] <- as.numeric(expr_dt[[cc]])
}

expr_dt <- expr_dt %>%
  mutate(
    WT_mean    = rowMeans(select(., `WT.1`, `WT.2`),             na.rm = TRUE),
    cpp8_mean  = rowMeans(select(., `TCX2-3-KO.1`,`TCX2-3-KO.2`), na.rm = TRUE),
    cpp11_mean = rowMeans(select(., `SOL1-8-KO.1`, `SOL1-8-KO.2`), na.rm = TRUE),
    log2FC_cpp8  = log2(cpp8_mean  + 1) - log2(WT_mean + 1),
    log2FC_cpp11 = log2(cpp11_mean + 1) - log2(WT_mean + 1)
  ) %>%
  filter(WT_mean > 0)   # 排除 WT 不表达基因，避免伪 FC

expr_fc <- expr_dt %>%
  select(chr, start, end, gene_id,
         WT_mean, cpp8_mean, cpp11_mean,
         log2FC_cpp8, log2FC_cpp11)

# =========================================================
# 8. TSS 映射到 1-kb repli bin
# =========================================================
expr_fc_dt <- as.data.table(expr_fc)
repli_dt   <- as.data.table(data_avg)

setkey(repli_dt,   chr, start, end)
setkey(expr_fc_dt, chr, start, end)

gene_bin_dt <- foverlaps(
  expr_fc_dt[, .(chr, start, end, gene_id,
                 WT_mean, cpp8_mean, cpp11_mean,
                 log2FC_cpp8, log2FC_cpp11)],
  repli_dt[,   .(chr, start, end,
                 WT_WRT, tcx2_3_WRT, sol1_8_WRT,
                 deltaRT_cpp8, deltaRT_cpp11)],
  nomatch = 0L
)

gene_bin_dt <- unique(gene_bin_dt[, .(
  gene_id, chr, start, end,
  WT_WRT, tcx2_3_WRT, sol1_8_WRT,
  deltaRT_cpp8, deltaRT_cpp11,
  WT_mean, cpp8_mean, cpp11_mean,
  log2FC_cpp8, log2FC_cpp11
)])

# =========================================================
# 9. 压缩为严格 gene-level（一基因一行）
# =========================================================
gene_bin_unique <- gene_bin_dt %>%
  as.data.frame() %>%
  group_by(gene_id) %>%
  summarise(
    deltaRT_cpp8  = mean(deltaRT_cpp8,  na.rm = TRUE),
    deltaRT_cpp11 = mean(deltaRT_cpp11, na.rm = TRUE),
    WT_WRT        = mean(WT_WRT,        na.rm = TRUE),
    tcx2_3_WRT    = mean(tcx2_3_WRT,    na.rm = TRUE),
    sol1_8_WRT    = mean(sol1_8_WRT,    na.rm = TRUE),
    WT_mean       = first(WT_mean),
    cpp8_mean     = first(cpp8_mean),
    cpp11_mean    = first(cpp11_mean),
    log2FC_cpp8   = first(log2FC_cpp8),
    log2FC_cpp11  = first(log2FC_cpp11),
    .groups = "drop"
  )

cat("gene-level行数:", nrow(gene_bin_unique),
    "| 唯一gene数:", length(unique(gene_bin_unique$gene_id)), "\n")

fwrite(as.data.table(gene_bin_unique),
       "gene_TSS_mapped_to_repli_bin.unique_gene_level.txt", sep = "\t")

# =========================================================
# 16. 模块A：按 WT WRT 五分位数分层分析
#     Q1 = 最早复制，Q5 = 最晚复制
# =========================================================

# --- 16.1 打 Q1-Q5 标签 ---
gene_bin_unique <- gene_bin_unique %>%
  mutate(
    RT_quintile = ntile(WT_WRT, 5),
    RT_quintile_label = factor(
      RT_quintile,
      levels = 1:5,
      labels = c("Q1 (Earliest)", "Q2", "Q3", "Q4", "Q5 (Latest)")
    )
  )

# --- 16.2 每个 quintile 内计算相关性 ---
quintile_cor <- gene_bin_unique %>%
  group_by(RT_quintile_label) %>%
  summarise(
    n                = n(),
    pearson_r_cpp8   = cor(deltaRT_cpp8,  log2FC_cpp8,  use = "complete.obs", method = "pearson"),
    pearson_p_cpp8   = tryCatch(cor.test(deltaRT_cpp8,  log2FC_cpp8,  method = "pearson")$p.value,  error = function(e) NA),
    spearman_r_cpp8  = cor(deltaRT_cpp8,  log2FC_cpp8,  use = "complete.obs", method = "spearman"),
    spearman_p_cpp8  = tryCatch(cor.test(deltaRT_cpp8,  log2FC_cpp8,  method = "spearman")$p.value, error = function(e) NA),
    pearson_r_cpp11  = cor(deltaRT_cpp11, log2FC_cpp11, use = "complete.obs", method = "pearson"),
    pearson_p_cpp11  = tryCatch(cor.test(deltaRT_cpp11, log2FC_cpp11, method = "pearson")$p.value,  error = function(e) NA),
    spearman_r_cpp11 = cor(deltaRT_cpp11, log2FC_cpp11, use = "complete.obs", method = "spearman"),
    spearman_p_cpp11 = tryCatch(cor.test(deltaRT_cpp11, log2FC_cpp11, method = "spearman")$p.value, error = function(e) NA),
    .groups = "drop"
  )

print(quintile_cor)
fwrite(as.data.table(quintile_cor), "RT_quintile_correlation_summary.txt", sep = "\t")

# --- 16.3 Facet 散点图（每个 quintile 一格）---
plot_facet_quintile <- function(df, xcol, ycol, mutant_label, out_pdf) {
  cor_labels <- df %>%
    group_by(RT_quintile_label) %>%
    summarise(
      n = n(),
      r = round(cor(.data[[xcol]], .data[[ycol]], use = "complete.obs", method = "pearson"), 3),
      p = signif(cor.test(.data[[xcol]], .data[[ycol]], method = "pearson")$p.value, 2),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("n=", n, "\nr=", r, "\nP=", p))
  
  p <- ggplot(df, aes_string(x = xcol, y = ycol)) +
    geom_point(alpha = 0.35, size = 0.8, color = "#2C7FB8") +
    geom_smooth(method = "lm", se = TRUE, color = "#D95F0E", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    facet_wrap(~ RT_quintile_label, nrow = 1, scales = "free") +
    geom_text(
      data = cor_labels,
      aes(label = label),
      x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
      size = 3, inherit.aes = FALSE
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "#EFF3FF"),
      strip.text       = element_text(face = "bold", size = 10),
      panel.spacing    = unit(0.8, "lines")
    ) +
    labs(
      x     = expression(Delta*RT~"(mutant - WT WRT)"),
      y     = "RNA log2FC",
      title = paste0(mutant_label, ": ΔRT vs RNA log2FC stratified by WT RT quintile")
    )
  
  ggsave(out_pdf, p, width = 16, height = 5, device = cairo_pdf)
  return(p)
}

p_quintile_cpp8 <- plot_facet_quintile(
  gene_bin_unique, "deltaRT_cpp8",  "log2FC_cpp8",
  "cpp8 (tcx2.3)", "cpp8_deltaRT_vs_log2FC_by_RT_quintile.pdf"
)
p_quintile_cpp11 <- plot_facet_quintile(
  gene_bin_unique, "deltaRT_cpp11", "log2FC_cpp11",
  "cpp11 (sol1.8)", "cpp11_deltaRT_vs_log2FC_by_RT_quintile.pdf"
)
print(p_quintile_cpp8)
print(p_quintile_cpp11)

# --- 16.4 相关性随 RT quintile 变化的趋势折线图 ---
cor_long <- quintile_cor %>%
  select(RT_quintile_label,
         cpp8 = pearson_r_cpp8, cpp11 = pearson_r_cpp11) %>%
  tidyr::pivot_longer(c(cpp8, cpp11), names_to = "mutant", values_to = "pearson_r") %>%
  mutate(mutant = recode(mutant,
                         "cpp8"  = "cpp8 (tcx2.3)",
                         "cpp11" = "cpp11 (sol1.8)"
  ))

p_cor_trend_quintile <- ggplot(cor_long,
                               aes(x = RT_quintile_label, y = pearson_r, color = mutant, group = mutant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("cpp8 (tcx2.3)" = "#D95F0E", "cpp11 (sol1.8)" = "#2C7FB8")) +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank()) +
  labs(
    x     = "WT RT Quintile (Q1=Earliest → Q5=Latest)",
    y     = "Pearson r (ΔRT vs RNA log2FC)",
    title = "Correlation between ΔRT and RNA log2FC across WT RT quintiles"
  )

ggsave("RT_quintile_correlation_trend.pdf",
       p_cor_trend_quintile, width = 8, height = 5, device = cairo_pdf)
print(p_cor_trend_quintile)

# =========================================================
# 17. 模块B：按 WT TPM 表达量分组分析
#     与文章 Fig.2D 保持一致的4档（T=0已在第7节排除）
# =========================================================

# --- 17.1 打表达量分组标签 ---
gene_bin_unique <- gene_bin_unique %>%
  mutate(
    expr_group = case_when(
      WT_mean > 0   & WT_mean <= 1   ~ "0<T≤1",
      WT_mean > 1   & WT_mean <= 10  ~ "1<T≤10",
      WT_mean > 10  & WT_mean <= 100 ~ "10<T≤100",
      WT_mean > 100                  ~ "T>100",
      TRUE                           ~ NA_character_
    ),
    expr_group = factor(expr_group,
                        levels = c("0<T≤1", "1<T≤10", "10<T≤100", "T>100"))
  )

print(table(gene_bin_unique$expr_group))

# --- 17.2 每个表达量组内计算相关性 ---
expr_group_cor <- gene_bin_unique %>%
  filter(!is.na(expr_group)) %>%
  group_by(expr_group) %>%
  summarise(
    n                = n(),
    pearson_r_cpp8   = cor(deltaRT_cpp8,  log2FC_cpp8,  use = "complete.obs", method = "pearson"),
    pearson_p_cpp8   = tryCatch(cor.test(deltaRT_cpp8,  log2FC_cpp8,  method = "pearson")$p.value,  error = function(e) NA),
    spearman_r_cpp8  = cor(deltaRT_cpp8,  log2FC_cpp8,  use = "complete.obs", method = "spearman"),
    spearman_p_cpp8  = tryCatch(cor.test(deltaRT_cpp8,  log2FC_cpp8,  method = "spearman")$p.value, error = function(e) NA),
    pearson_r_cpp11  = cor(deltaRT_cpp11, log2FC_cpp11, use = "complete.obs", method = "pearson"),
    pearson_p_cpp11  = tryCatch(cor.test(deltaRT_cpp11, log2FC_cpp11, method = "pearson")$p.value,  error = function(e) NA),
    spearman_r_cpp11 = cor(deltaRT_cpp11, log2FC_cpp11, use = "complete.obs", method = "spearman"),
    spearman_p_cpp11 = tryCatch(cor.test(deltaRT_cpp11, log2FC_cpp11, method = "spearman")$p.value, error = function(e) NA),
    .groups = "drop"
  )

print(expr_group_cor)
fwrite(as.data.table(expr_group_cor), "expr_group_correlation_summary.txt", sep = "\t")

# --- 17.3 Facet 散点图（每个表达量组一格）---
plot_facet_expr <- function(df, xcol, ycol, mutant_label, out_pdf) {
  df_plot <- df %>% filter(!is.na(expr_group))
  
  cor_labels <- df_plot %>%
    group_by(expr_group) %>%
    summarise(
      n = n(),
      r = round(cor(.data[[xcol]], .data[[ycol]], use = "complete.obs", method = "pearson"), 3),
      p = signif(tryCatch(
        cor.test(.data[[xcol]], .data[[ycol]], method = "pearson")$p.value,
        error = function(e) NA), 2),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("n=", n, "\nr=", r, "\nP=", p))
  
  p <- ggplot(df_plot, aes_string(x = xcol, y = ycol)) +
    geom_point(alpha = 0.35, size = 0.8, color = "#2C7FB8") +
    geom_smooth(method = "lm", se = TRUE, color = "#D95F0E", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    facet_wrap(~ expr_group, nrow = 1, scales = "free") +
    geom_text(
      data = cor_labels,
      aes(label = label),
      x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
      size = 3, inherit.aes = FALSE
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "#F0F7E6"),
      strip.text       = element_text(face = "bold", size = 10),
      panel.spacing    = unit(0.8, "lines")
    ) +
    labs(
      x     = expression(Delta*RT~"(mutant - WT WRT)"),
      y     = "RNA log2FC",
      title = paste0(mutant_label, ": ΔRT vs RNA log2FC stratified by WT expression level")
    )
  
  ggsave(out_pdf, p, width = 14, height = 5, device = cairo_pdf)
  return(p)
}

p_expr_cpp8 <- plot_facet_expr(
  gene_bin_unique, "deltaRT_cpp8",  "log2FC_cpp8",
  "cpp8 (tcx2.3)", "cpp8_deltaRT_vs_log2FC_by_expr_group.pdf"
)
p_expr_cpp11 <- plot_facet_expr(
  gene_bin_unique, "deltaRT_cpp11", "log2FC_cpp11",
  "cpp11 (sol1.8)", "cpp11_deltaRT_vs_log2FC_by_expr_group.pdf"
)
print(p_expr_cpp8)
print(p_expr_cpp11)

# --- 17.4 相关性随表达量组变化的趋势折线图 ---
expr_cor_long <- expr_group_cor %>%
  select(expr_group,
         cpp8 = pearson_r_cpp8, cpp11 = pearson_r_cpp11) %>%
  tidyr::pivot_longer(c(cpp8, cpp11), names_to = "mutant", values_to = "pearson_r") %>%
  mutate(mutant = recode(mutant,
                         "cpp8"  = "cpp8 (tcx2.3)",
                         "cpp11" = "cpp11 (sol1.8)"
  ))

p_cor_trend_expr <- ggplot(expr_cor_long,
                           aes(x = expr_group, y = pearson_r, color = mutant, group = mutant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("cpp8 (tcx2.3)" = "#D95F0E", "cpp11 (sol1.8)" = "#2C7FB8")) +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank()) +
  labs(
    x     = "WT Expression Level (TPM)",
    y     = "Pearson r (ΔRT vs RNA log2FC)",
    title = "Correlation between ΔRT and RNA log2FC across expression groups"
  )

ggsave("expr_group_correlation_trend.pdf",
       p_cor_trend_expr, width = 8, height = 5, device = cairo_pdf)
print(p_cor_trend_expr)



