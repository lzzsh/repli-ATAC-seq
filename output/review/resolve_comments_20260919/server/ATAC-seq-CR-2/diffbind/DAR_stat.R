# ================== 加载包 ==================
library(dplyr)
library(ggplot2)
library(ChIPseeker)
library(GenomicFeatures)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/")

# ================== 读入所有样本（文件名不变） ==================
files <- list(
  CPP8_1  = "TCX2_1_vs_WT_edgeR_sig.bed",
  CPP8_3  = "TCX2_3_vs_WT_edgeR_sig.bed",
  CPP11_5 = "sol1_5_vs_WT_edgeR_sig.bed",
  CPP11_8 = "sol1_8_vs_WT_edgeR_sig.bed"
)

# ================== Step1: Up / Down 统计 ==================
updown_stats <- lapply(names(files), function(nm){
  df <- read.table(files[[nm]], header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "strand", "log2FC")
  
  if(!"log2FC" %in% colnames(df)) stop(paste("列 log2FC 在文件中不存在:", nm))
  
  df %>%
    mutate(Regulation = ifelse(log2FC > 0, "Up", "Down")) %>%
    group_by(Regulation) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(Sample = nm)
}) %>% bind_rows()

# ---- 绘制 Up/Down 柱状图（Sample 标签斜体）----
p1 <- ggplot(updown_stats, aes(x = Sample, y = n, fill = Regulation)) +
  geom_col(position = "dodge", color = "black", width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Up"="#E69191", "Down"="#92B5CA")) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "italic")  # 斜体
  ) +
  labs(title = "Up/Down Peaks per condition", x = NULL, y = "Peak Count")

ggsave("UpDown_peak_count.pdf", p1, width = 6, height = 4.5)


# ================== Step2: 基因组注释分布 ==================
txdb <- makeTxDbFromGFF(
  "/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3"
)

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
    summarise(n = n(), .groups = "drop") %>%
    mutate(Sample = nm)
}) %>% bind_rows()

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

genomic_dist_list$annotation_group <- factor(genomic_dist_list$annotation_group,
                                             levels = feature_order)

my_colors <- c("#C25759", "#E69191", "#599CB4", "#92B5CA")

# ---- 柱状图：legend + x轴样本标签斜体 ----
p2_bar <- ggplot(genomic_dist_list,
                 aes(x = annotation_group, y = n, fill = Sample)) +
  geom_col(position = "dodge", color = "black") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(face = "italic"),   # 斜体
    legend.title = element_blank()
  ) +
  labs(title = "Genomic Distribution of Differential Peaks",
       x = "", y = "Count") +
  scale_fill_manual(values = my_colors)

ggsave("Genomic_distribution_all_samples_promoter_subtype.pdf",
       p2_bar, width = 6, height = 4.5)

# ---- 折线图：legend + x轴样本标签斜体 ----
p2_line <- ggplot(genomic_dist_list,
                  aes(x = annotation_group, y = n, group = Sample, color = Sample)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = my_colors) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(face = "italic"),   # 斜体
    legend.title = element_blank()
  ) +
  labs(title = "Genomic Distribution of Differential Peaks",
       x = "", y = "Peak Count")

p2_line


# ================== Venn：把 TCX2/SOL1 名称换成 CPP8/CPP11 + 标签斜体 ==================
library(VennDiagram)

peak_list <- lapply(files, function(f){
  df <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "strand", "log2FC")
  df$peak_id <- paste(df$chr, df$start, df$end, sep = "_")
  df$peak_id
})
names(peak_list) <- names(files)

sapply(peak_list, length)

shared_all   <- Reduce(intersect, peak_list)
shared_CPP8  <- intersect(peak_list$CPP8_1,  peak_list$CPP8_3)
shared_CPP11 <- intersect(peak_list$CPP11_5, peak_list$CPP11_8)

cat("CPP8_1 ∩ CPP8_3 =", length(shared_CPP8), "\n")
cat("CPP11_5 ∩ CPP11_8 =", length(shared_CPP11), "\n")
cat("所有四组共有 =", length(shared_all), "\n")

venn.plot <- venn.diagram(
  x = list(
    CPP8_1  = peak_list$CPP8_1,
    CPP8_3  = peak_list$CPP8_3,
    CPP11_5 = peak_list$CPP11_5,
    CPP11_8 = peak_list$CPP11_8
  ),
  filename = "Venn_Differential_Peaks.tiff",
  
  main = "Overlap of Differential Peaks",
  main.cex = 1.5,
  main.fontface = "bold.italic",   # 标题：斜体 + 粗体
  main.fontfamily = "Arial",
  
  # 把 _ 改成 -
  category.names = gsub("_", "-", names(files)),
  
  col = "black",
  fill = c("#C25759", "#E69191", "#599CB4", "#92B5CA"),
  alpha = 0.5,
  
  cex = 1.2,
  fontface = "plain",              # 中间数字：普通体
  fontfamily = "Arial",
  
  cat.cex = 1.2,
  cat.fontface = "bold.italic",    # 标签：斜体 + 粗体
  cat.fontfamily = "Arial",
  
  cat.pos = 0,
  margin = 0.05,
  lwd = 2
)
