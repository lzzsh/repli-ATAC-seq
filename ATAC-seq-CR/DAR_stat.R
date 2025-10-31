# ================== 加载包 ==================
library(dplyr)
library(ggplot2)
library(ChIPseeker)
library(GenomicFeatures)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/")

# ================== 读入所有样本 ==================
files <- list(
  TCX2_1 = "TCX2_1_vs_WT_edgeR_sig.bed",
  TCX2_3 = "TCX2_3_vs_WT_edgeR_sig.bed",
  SOL1_5 = "sol1_5_vs_WT_edgeR_sig.bed",
  SOL1_8 = "sol1_8_vs_WT_edgeR_sig.bed"
)

# ================== Step1: Up / Down 统计 ==================
updown_stats <- lapply(names(files), function(nm){
  df <- read.table(files[[nm]], header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "strand", "log2FC")
  
  # 检查 log2FC 列是否存在
  if(!"log2FC" %in% colnames(df)) {
    stop(paste("列 log2FC 在文件中不存在:", nm))
  }
  
  df %>%
    mutate(Regulation = ifelse(log2FC > 0, "Up", "Down")) %>%
    group_by(Regulation) %>%
    summarise(n = n()) %>%
    mutate(Sample = nm)
}) %>%
  bind_rows()

# ---- 绘制 Up/Down 柱状图 ----
p1 <- ggplot(updown_stats, aes(x = Sample, y = n, fill = Regulation)) +
  geom_col(position = "dodge", color = "black", width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7), 
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Up"="#E69191", "Down"="#92B5CA")) +
  theme_classic(base_size = 14) +
  labs(title = "Up/Down Peaks per condition", x = NULL, y = "Peak Count")

ggsave("UpDown_peak_count.pdf", p1, width = 6, height = 4.5)


# ================== Step2: 基因组注释分布 ==================
txdb <- makeTxDbFromGFF(
  "/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3"
)

# ===================== 基因组注释 + Promoter细分 =====================
genomic_dist_list <- lapply(names(files), function(nm){
  df <- read.table(files[[nm]], header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "strand", "log2FC")
  
  gr <- GRanges(seqnames = df$chr,
                ranges = IRanges(start = df$start, end = df$end),
                strand = df$strand)
  
  peakAnno <- annotatePeak(gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
  anno_df <- as.data.frame(peakAnno)
  
  anno_df %>%
    mutate(annotation_group = case_when(
      # promoter 区域进一步细分
      grepl("Promoter", annotation, ignore.case = TRUE) & distanceToTSS <= 1000  ~ "Promoter (≤1kb)",
      grepl("Promoter", annotation, ignore.case = TRUE) & distanceToTSS > 1000 & distanceToTSS <= 2000 ~ "Promoter (1–2kb)",
      grepl("Promoter", annotation, ignore.case = TRUE) & distanceToTSS > 2000 & distanceToTSS <= 3000 ~ "Promoter (2–3kb)",
      grepl("Promoter", annotation, ignore.case = TRUE) ~ "Promoter (>3kb)",
      grepl("Intron", annotation, ignore.case = TRUE) ~ "Intron",
      grepl("Exon", annotation, ignore.case = TRUE) ~ "Exon",
      grepl("3'UTR", annotation, ignore.case = TRUE) ~ "3'UTR",
      grepl("5'UTR", annotation, ignore.case = TRUE) ~ "5'UTR",
      grepl("Intergenic", annotation, ignore.case = TRUE) ~ "Intergenic",
      TRUE ~ "Other"
    )) %>%
    group_by(annotation_group) %>%
    summarise(n = n()) %>%
    mutate(Sample = nm)
}) %>%
  bind_rows()

feature_order <- c("Promoter (≤1kb)", 
                   "Promoter (1–2kb)", 
                   "Promoter (2–3kb)", 
                   "Promoter (>3kb)",
                   "5'UTR", 
                   "Exon", 
                   "Intron", 
                   "3'UTR", 
                   "Intergenic", 
                   "Other")

# 将 annotation_group 转换为 factor
genomic_dist_list$annotation_group <- factor(genomic_dist_list$annotation_group,
                                             levels = feature_order)

# ===================== 绘图 =====================
my_colors <- c("#C25759", "#E69191", "#599CB4", "#92B5CA")

p2 <- ggplot(genomic_dist_list,
             aes(x = annotation_group, y = n, fill = Sample)) +
  geom_col(position = "dodge", color = "black") +
  theme_classic(base_size = 13) +
  labs(title = "Genomic Distribution of Differential Peaks",
       x = "", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = my_colors)

ggsave("Genomic_distribution_all_samples_promoter_subtype.pdf", p2, width = 6, height = 4.5)

p2 <- ggplot(genomic_dist_list,
             aes(x = annotation_group, y = n, group = Sample, color = Sample)) +
  geom_line(size = 1) +                                    # 连线
  geom_point(size = 3) +                                   # 点
  scale_color_manual(values = my_colors) +                 # 手动配色
  theme_classic(base_size = 13) +
  labs(title = "Genomic Distribution of Differential Peaks",
       x = "", y = "Peak Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())
p2


############
# ================== 加载包 ==================
library(dplyr)
library(ggplot2)
library(VennDiagram)  # ← 关键包

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/")

# ================== 读入所有样本 ==================
files <- list(
  TCX2_1 = "TCX2_1_vs_WT_edgeR_sig.bed",
  TCX2_3 = "TCX2_3_vs_WT_edgeR_sig.bed",
  SOL1_5 = "sol1_5_vs_WT_edgeR_sig.bed",
  SOL1_8 = "sol1_8_vs_WT_edgeR_sig.bed"
)

# ---- 1. 读取所有 peaks 并生成唯一ID ----
peak_list <- lapply(files, function(f){
  df <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "strand", "log2FC")
  df$peak_id <- paste(df$chr, df$start, df$end, sep = "_")  # 唯一ID
  df$peak_id
})
names(peak_list) <- names(files)

# ---- 2. 查看每个样本的峰数量 ----
sapply(peak_list, length)

# ---- 3. 计算交集 ----
shared_all  <- Reduce(intersect, peak_list)
shared_TCX2 <- intersect(peak_list$TCX2_1, peak_list$TCX2_3)
shared_SOL1 <- intersect(peak_list$SOL1_5, peak_list$SOL1_8)

cat("TCX2_1 ∩ TCX2_3 =", length(shared_TCX2), "\n")
cat("SOL1_5 ∩ SOL1_8 =", length(shared_SOL1), "\n")
cat("所有四组共有 =", length(shared_all), "\n")

# ---- 4. 绘制四集合 Venn 图 ----
venn.plot <- venn.diagram(
  x = list(
    TCX2_1 = peak_list$TCX2_1,
    TCX2_3 = peak_list$TCX2_3,
    SOL1_5 = peak_list$SOL1_5,
    SOL1_8 = peak_list$SOL1_8
  ),
  filename = "Venn_Differential_Peaks.tiff",  # 输出文件，可改为NULL直接显示
  main = "Overlap of Differential Peaks",
  main.cex = 1.5,
  category.names = names(files),
  col = "black",                      # 圆圈边框
  fill = c("#C25759", "#E69191", "#599CB4", "#92B5CA"),  # 填充色
  alpha = 0.5,                        # 透明度
  cex = 1.2,                          # 数字字体大小
  cat.cex = 1.2,                      # 标签字体大小
  cat.fontface = "bold",              # 标签加粗
  cat.pos = 0,                        # 自动调整标签位置
  margin = 0.05,
  lwd = 2                             # 圆圈线宽
)

