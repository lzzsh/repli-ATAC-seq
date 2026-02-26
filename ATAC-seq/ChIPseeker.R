library(ChIPseeker)
library(GenomicFeatures)
library(readxl)

setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR-2/diffbind/")

# Build TxDb from rice GFF3 annotation
ricegff <- makeTxDbFromGFF("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3")

# annotate peaks
gene_anno <- read_excel("../../ATAC-seq-CR/diffbind/gene_function.xlsx")
colnames(gene_anno)[1] <- "geneId" 
deseq2_peaks <- read.table("./TCX2_3_vs_WT_deseq2_sig.bed")
edgeR_peaks <- read.table("./TCX2_3_vs_WT_edgeR_sig.bed")

# 给列命名（如果你的 bed 文件确实是这 5 列）
colnames(deseq2_peaks) <- c("chr","start","end","strand","log2FC")
colnames(edgeR_peaks) <- c("chr","start","end","strand","log2FC")

# 2. 转换为 GRanges
gr_deseq2 <- GRanges(seqnames = edgeR_peaks$chr,
                     ranges   = IRanges(start = edgeR_peaks$start,
                                        end   = edgeR_peaks$end),
                     strand   = edgeR_peaks$strand,
                     log2FC   = edgeR_peaks$log2FC)

# tcx2 binding sites
tcx2_cuttag <- read.table("/storage2/liuxiaodongLab/liaozizhuo/Projects/cuttag_tcx/macs2/macs2_p1e-5/tcx_cut_out/tcx_cut_peaks_sort.narrowPeak")

colnames(tcx2_cuttag) <- c(
  "chr", "start", "end", "name", "score", "strand",
  "signalValue", "pValue", "qValue", "peak"
)

gr_cuttag <- GRanges(
  seqnames = tcx2_cuttag$chr,
  ranges   = IRanges(start = tcx2_cuttag$start, end = tcx2_cuttag$end)
)

gr_deseq2_expanded <- promoters(gr_deseq2, upstream = 10000, downstream = width(gr_deseq2) + 10000)

overlap_hits <- findOverlaps(gr_deseq2_expanded, gr_cuttag)
gr_deseq2_with_cuttag <- gr_deseq2[unique(queryHits(overlap_hits))]

peakAnno <- annotatePeak(
  gr_deseq2_with_cuttag,
  TxDb = ricegff,
  tssRegion = c(-3000, 3000),
  annoDb = NULL
)

# 转为 data.frame
peakAnno_df <- as.data.frame(peakAnno)
peakAnno_df <- merge(
  peakAnno_df,
  gene_anno[,c(1,8)],
  by = "geneId",
  all.x = TRUE
)

# 筛选距离 TSS 小于 3kb 的 peaks
peakAnno_nearTSS <- subset(peakAnno_df, abs(distanceToTSS) < 3000)

peakAnno_nearTSS_sorted <- peakAnno_nearTSS[order(-abs(peakAnno_nearTSS$log2FC)), ]

# 查看前几行
head(peakAnno_nearTSS_sorted)

write.csv(peakAnno_nearTSS_sorted,
            file = "./TCX2-3_CR_WT_ES_edger_sig_with_cuttag_nearTSS3kb_sorted.csv",
            quote = FALSE, col.names = TRUE, row.names = FALSE)


############################################################
# 新增：检验 TCX2 CUT&Tag 在 DAR 中的富集情况
############################################################

library(GenomicRanges)
library(dplyr)

# gr_deseq2 = 差异ATAC peaks (edgeR sig)
# gr_cuttag = TCX2 CUT&Tag peaks

# -------------------------
# 1) 构建所有DAR是否被TCX2结合的标记
# -------------------------
overlap_all <- findOverlaps(gr_deseq2, gr_cuttag)
gr_deseq2$TCX2_bound <- seq_along(gr_deseq2) %in% queryHits(overlap_all)

# -------------------------
# 2) 构建背景 = 所有 ATAC peaks（需要你全peak集）
# 如果你有全peak文件，把它改成你自己的路径
# -------------------------
all_peaks <- read.table("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT_all_org.gff3")
colnames(all_peaks) <- c("chr","start","end","RT")

gr_all <- GRanges(
  seqnames = all_peaks$chr,
  ranges   = IRanges(start = all_peaks$start, end = all_peaks$end)
)

# 标记所有peak是否为DAR
hits_dar <- findOverlaps(gr_all, gr_deseq2)
gr_all$isDAR <- seq_along(gr_all) %in% queryHits(hits_dar)

# 标记所有peak是否TCX2结合
hits_cuttag <- findOverlaps(gr_all, gr_cuttag)
gr_all$TCX2_bound <- seq_along(gr_all) %in% queryHits(hits_cuttag)

# -------------------------
# 3) 构建 2x2 表格: Bound vs DAR
# -------------------------
tbl <- table(gr_all$TCX2_bound, gr_all$isDAR)
print(tbl)
#           Non-DAR   DAR
# Unbound      d       c
# Bound        b       a

# -------------------------
# 4) Fisher Test
# -------------------------
res <- fisher.test(tbl)
print(res)
