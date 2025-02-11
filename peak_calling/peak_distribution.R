# 加载必要的 R 包
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

# 1 读取 GFF3 文件，获取基因注释信息
gff_file <- "/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3"
gff_data <- import(gff_file, format="gff3")

# 2 读取 WOX11 CUT&Tag peak 数据
peaks_file <- "/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_cut_peaks_top1000.narrowPeak"
peaks <- read.table(peaks_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(peaks) <- c("chr", "start", "end", "V4", "score", "strand", "signalValue", "pValue", "qValue", "peak")

# 转换 peaks 为 GenomicRanges
gr_peaks <- GRanges(seqnames = peaks$chr, 
                    ranges = IRanges(start = peaks$start, end = peaks$end))

# 3 提取不同基因区域
# 提取基因 (gene)、mRNA (mRNA)、5'UTR、3'UTR
genes <- gff_data[gff_data$type == "gene"]
mrnas <- gff_data[gff_data$type == "mRNA"]
utr5 <- gff_data[gff_data$type == "five_prime_UTR"]
utr3 <- gff_data[gff_data$type == "three_prime_UTR"]

# 4 计算启动子区域（基因起始点上游 2000 bp）
promoters <- genes
start(promoters) <- start(promoters) - 2000  # 启动子上游 2000 bp
end(promoters) <- start(promoters) + 2000   # 启动子范围 4000 bp

# 5 计算 peaks 在不同区域的交集
overlap_promoter <- countOverlaps(gr_peaks, promoters) > 0
overlap_utr5 <- countOverlaps(gr_peaks, utr5) > 0
overlap_utr3 <- countOverlaps(gr_peaks, utr3) > 0
overlap_gene <- countOverlaps(gr_peaks, genes) > 0

# 计算基因间区（没有匹配到基因的 peaks 视为基因间区）
overlap_intergenic <- !(overlap_promoter | overlap_utr5 | overlap_utr3 | overlap_gene)

# 6 统计 peaks 在各个区域的数量
peak_counts <- data.frame(
  Category = c("Promoter", "5' UTR", "3' UTR", "Gene", "Intergenic"),
  Count = c(sum(overlap_promoter), sum(overlap_utr5), sum(overlap_utr3), sum(overlap_gene), sum(overlap_intergenic))
)

# 计算百分比
peak_counts <- peak_counts %>%
  mutate(Percentage = Count / sum(Count) * 100)

# 7 显示结果
print(peak_counts)

# 8 绘制直方图
ggplot(peak_counts, aes(x=Category, y=Percentage, fill=Category)) +
  geom_bar(stat="identity", width = 0.9, color = "black") +
  scale_fill_manual(values=c("#D4E4C3", "#EAC9C1", "#F0DBCD", "#D6C4C4", "#E8E6C4")) +
  theme_minimal() +
  labs(title="WOX11 CUT&Tag Peaks Distribution", y="Percentage (%)", x="Genomic Region") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme_classic()
