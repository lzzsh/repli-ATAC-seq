# ================== 加载包 ==================
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(ggplot2)
library(dplyr)

# ================== Step1: 构建 TxDb (用水稻 gff3) ==================
txdb <- makeTxDbFromGFF("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3")

# ================== Step2: 读入差异 Peaks (以 DESeq2 结果为例) ==================
# bed 文件包含: chr, start, end, strand, log2FC
deseq2_peaks <- read.table("./WT_CR_ES_deseq2_sig.bed")
colnames(deseq2_peaks) <- c("chr","start","end","strand","log2FC")

# 转成 GRanges
gr_deseq2 <- GRanges(
  seqnames = deseq2_peaks$chr,
  ranges   = IRanges(start = deseq2_peaks$start, end = deseq2_peaks$end),
  strand   = deseq2_peaks$strand,
  log2FC   = deseq2_peaks$log2FC
)

# ================== Step3: 注释 Peaks ==================
peakAnno <- annotatePeak(
  gr_deseq2,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = NULL
)

# 转成 data.frame
peakAnno_df <- as.data.frame(peakAnno)

# ================== Step4: Fig4B (Up 在上，Down 在下) ==================
peakAnno_df <- peakAnno_df %>%
  mutate(Regulation = ifelse(log2FC > 0, "Up", "Down"))

# 统计每类的数量
count_df <- peakAnno_df %>%
  group_by(Regulation) %>%
  summarise(n = n())

p <- ggplot(count_df, aes(x = Regulation, y = n, fill = Regulation)) +
  geom_col(width = 0.8, color="black") +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_manual(values = c("Up"="firebrick","Down"="steelblue")) +
  labs(title="Differential Peaks (DESeq2)", x=NULL, y="Count") +
  theme_classic() +
  theme(legend.position="none")
ggsave("./DP_tcx2-1_baplot.pdf", p, width = 4, height = 6)

# ================== Step5: Fig4C (基因组注释分布柱状图) ==================
anno_count <- peakAnno_df %>%
  filter(!grepl("Downstream", annotation, ignore.case = TRUE)) %>%
  mutate(annotation_group = case_when(
    grepl("Promoter", annotation, ignore.case = TRUE) ~ annotation,  # 保留 promoter 分段
    grepl("Intron", annotation, ignore.case = TRUE) ~ "Intron",
    grepl("Exon", annotation, ignore.case = TRUE) ~ "Exon",
    grepl("3'UTR", annotation, ignore.case = TRUE) ~ "3'UTR",
    grepl("5'UTR", annotation, ignore.case = TRUE) ~ "5'UTR",
    grepl("Intergenic", annotation, ignore.case = TRUE) ~ "Intergenic",
    TRUE ~ annotation
  )) %>%
  group_by(annotation_group) %>%
  summarise(count = n())

ggplot(anno_count, aes(x = reorder(annotation_group, -count), y = count)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_bw() +
  labs(title="Genomic Distribution of Differential Peaks",
       x="Genomic Feature", y="Count") +
  theme_classic()
