library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/")

atac_deg <- read.csv("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/TCX2-3_CR_WT_ES_edger_sig_with_cuttag_nearTSS3kb_sorted.csv")
rna_deg <- read.csv("./edgeR_TCX2-3-KO_vs_WT.csv")

colnames(rna_deg)[1] <- "MSU"
colnames(atac_deg)[1] <- "MSU"

atac_deg_common <- atac_deg %>%
  inner_join(rna_deg %>% dplyr::select(MSU, logFC), by = "MSU")
names(atac_deg_common)[names(atac_deg_common) == "log2FC"] <- "log2FC_ATAC"
names(atac_deg_common)[names(atac_deg_common) == "logFC"]  <- "log2FC_RNA"

atac_deg_common <- atac_deg_common %>%
  mutate(
    Regulation = case_when(
      log2FC_ATAC > 0 & log2FC_RNA > 0 ~ "Both_Up",
      log2FC_ATAC < 0 & log2FC_RNA < 0 ~ "Both_Down",
      TRUE ~ "Opposite"
    ),
    Regulation = factor(Regulation, levels = c("Both_Up", "Both_Down", "Opposite"))
  )

fit <- lm(log2FC_ATAC ~ log2FC_RNA, data = atac_deg_common)
r2 <- round(summary(fit)$r.squared, 2)

xlim_val <- ceiling(max(abs(atac_deg_common$log2FC_RNA), na.rm = TRUE))
ylim_val <- ceiling(max(abs(atac_deg_common$log2FC_ATAC), na.rm = TRUE))

p <- ggplot(atac_deg_common, aes(x = log2FC_RNA, y = log2FC_ATAC, color = Regulation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "#C25759", fill = "#E69191", size = 0.8, alpha = 0.25, fullrange = TRUE) +
  annotate("text", x = -xlim_val + 0.3, y = ylim_val - 0.3, label = sprintf("R² = %.2f", r2),
           hjust = 0, vjust = 1, size = 4.8, color = "black", fontface = "bold") +
  scale_color_manual(values = c("Both_Up"="#E69191","Both_Down"="#599CB4","Opposite"="grey70")) +
  labs(x = "RNA log2(Fold Change)", y = "ATAC log2(Fold Change)", color = "Regulation") +
  scale_x_continuous(limits = c(-xlim_val, xlim_val), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-ylim_val, ylim_val), expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.line = element_line(size = 0.6, color = "black"))

ggsave("Fig4F_ATAC_vs_RNA_with_fit.pdf", p, width = 6, height = 6)



########
# atac_deg_common <- atac_deg %>%
#   inner_join(rna_deg %>% dplyr::select(MSU, logFC), by = "MSU") 
# names(atac_deg_common)[names(atac_deg_common) == "log2FC"] <- "log2FC_ATAC"
# names(atac_deg_common)[names(atac_deg_common) == "logFC"]  <- "log2FC_RNA"
# 
# atac_deg_common <- atac_deg_common %>%
#   mutate(
#     RegulationType = case_when(
#       (log2FC_ATAC > 0 & log2FC_RNA > 0) | (log2FC_ATAC < 0 & log2FC_RNA < 0) ~ "Direct",
#       (log2FC_ATAC > 0 & log2FC_RNA < 0) | (log2FC_ATAC < 0 & log2FC_RNA > 0) ~ "Indirect",
#       TRUE ~ NA_character_
#     )
#   )
# 
# write.csv(atac_deg_common, "./ATAC_RNA_common_genes_TCX2-3.csv", row.names = FALSE)

library(VennDiagram)
library(dplyr)

# 读取数据
atac_deg <- read.csv("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/TCX2-3_CR_WT_ES_edger_sig_with_cuttag_nearTSS3kb_sorted.csv")
rna_deg  <- read.csv("./edgeR_TCX2-3-KO_vs_WT.csv")

colnames(atac_deg)[1] <- "MSU"
colnames(rna_deg)[1]  <- "MSU"

# 按 fold change 筛选
atac_deg_filt <- atac_deg %>% filter(abs(log2FC) > 0.3)
rna_deg_filt  <- rna_deg  %>% filter(abs(logFC)  > 1)

# 合并得到重叠基因信息
atac_rna_merge <- inner_join(
  atac_deg_filt %>% dplyr::select(MSU, log2FC_ATAC = log2FC),
  rna_deg_filt  %>% dplyr::select(MSU, log2FC_RNA = logFC),
  by = "MSU"
)

# 标注 direct / indirect
atac_rna_merge <- atac_rna_merge %>%
  mutate(
    RegulationType = case_when(
      (log2FC_ATAC > 0 & log2FC_RNA > 0) | (log2FC_ATAC < 0 & log2FC_RNA < 0) ~ "Direct",
      (log2FC_ATAC > 0 & log2FC_RNA < 0) | (log2FC_ATAC < 0 & log2FC_RNA > 0) ~ "Indirect",
      TRUE ~ NA_character_
    )
  )

# 获取基因集合
atac_genes <- unique(atac_deg_filt$MSU)
rna_genes  <- unique(rna_deg_filt$MSU)
common_genes <- unique(atac_rna_merge$MSU)

# 绘制Venn图
venn.plot <- venn.diagram(
  x = list(ATAC_DEG = atac_genes, RNA_DEG = rna_genes),
  filename = "Venn_ATAC_RNA_DEG.tiff",
  main = "Overlap between ATAC and RNA DEGs (filtered)",
  main.cex = 1.5,
  col = "black",
  fill = c("#599CB4", "#C25759"),
  alpha = 0.5,
  lwd = 2,
  cex = 1.2,
  cat.cex = 1.3,
  cat.fontface = "bold",
  cat.pos = c(0, 0),
  cat.dist = 0.02,
  margin = 0.1
)

# 加入注释
gene_anno <- read.csv("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/gene_annot.csv")
gene_anno$Msu <- trimws(gene_anno$Msu)

atac_rna_merge <- atac_rna_merge %>%
  distinct(MSU, .keep_all = TRUE)

common_anno <- gene_anno %>%
  filter(Msu %in% common_genes) %>%
  left_join(atac_rna_merge, by = c("Msu" = "MSU"))


# 导出
write.table(
  common_anno,
  "/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/common_genes_with_regulation.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)


library(dplyr)
library(stringr)

# 1. 读取数据 ------------------------------------------------
atac_deg <- read.csv("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/TCX2-3_CR_WT_ES_edger_sig_with_cuttag_nearTSS3kb_sorted.csv")
rna_deg  <- read.csv("./edgeR_TCX2-3-KO_vs_WT.csv")
gene_anno <- read.csv("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/gene_annot.csv")

# 统一列名
colnames(atac_deg)[1] <- "MSU"
colnames(rna_deg)[1]  <- "MSU"
gene_anno$Msu <- trimws(gene_anno$Msu)

# 2. 筛选显著差异基因 ----------------------------------------
atac_deg_filt <- atac_deg %>% filter(abs(log2FC) > 0.3)
rna_deg_filt  <- rna_deg  %>% filter(abs(logFC)  > 1)

# 3. 合并 ATAC 与 RNA DEG ------------------------------------
atac_rna_merge <- inner_join(
  atac_deg_filt %>% dplyr::select(MSU, log2FC_ATAC = log2FC),
  rna_deg_filt  %>% dplyr::select(MSU, log2FC_RNA = logFC),
  by = "MSU"
) %>%
  mutate(
    RegulationType = case_when(
      (log2FC_ATAC > 0 & log2FC_RNA > 0) | (log2FC_ATAC < 0 & log2FC_RNA < 0) ~ "Direct",
      (log2FC_ATAC > 0 & log2FC_RNA < 0) | (log2FC_ATAC < 0 & log2FC_RNA > 0) ~ "Indirect",
      TRUE ~ NA_character_
    )
  ) %>%
  distinct(MSU, .keep_all = TRUE)

# 4. 清理注释信息 --------------------------------------------
gene_anno_clean <- gene_anno %>%
  mutate(across(
    where(is.character),
    ~ str_replace_all(., "[^[:alnum:] ,.;:_-]", "")
  ))

# 5. 合并注释与调控信息 --------------------------------------
common_anno <- gene_anno_clean %>%
  filter(Msu %in% atac_rna_merge$MSU) %>%
  left_join(atac_rna_merge, by = c("Msu" = "MSU")) %>%
  mutate(
    log2FC_RNA = round(log2FC_RNA, 3),
    log2FC_ATAC = round(log2FC_ATAC, 3),
    # 若 Qsymbols 缺失或为空，则使用 Msu
    DisplayName = ifelse(!is.na(Qsymbols) & Qsymbols != "", Qsymbols, Msu),
    # 节点颜色映射
    NodeColor = case_when(
      RegulationType == "Direct" ~ "#C25759",   # 红色
      RegulationType == "Indirect" ~ "#599CB4", # 蓝色
      TRUE ~ "#D3D3D3"                          # 灰色
    )
  ) %>%
  filter(!is.na(Qsymbols) & Qsymbols != "NA")

# 6. 生成 Cytoscape Edge 文件 --------------------------------
edges <- common_anno %>%
  mutate(
    Source = "TCX2",
    Target = DisplayName,
    InteractionType = RegulationType
  ) %>%
  dplyr::select(Source, Target, InteractionType) %>%
  drop_na(Target)  # 删除 Target 为 NA 的行

write.table(
  edges,
  "/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/cytoscape_edges.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 7. 生成 Cytoscape Node 文件 --------------------------------
nodes <- common_anno %>%
  mutate(ID = DisplayName, NodeType = "Gene") %>%
  dplyr::select(ID, log2FC_RNA, log2FC_ATAC, RegulationType, NodeType, NodeColor, everything()) %>%
  drop_na(ID)  # 删除 ID 为 NA 的行

write.table(
  nodes,
  "/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/cytoscape_nodes.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)
