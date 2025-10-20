library(ChIPseeker)
library(GenomicFeatures)
library(readxl)

# Build TxDb from rice GFF3 annotation
ricegff <- makeTxDbFromGFF("/storage2/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3")

# annotate peaks
gene_anno <- read_excel("./gene_function.xlsx")
colnames(gene_anno)[1] <- "geneId" 
deseq2_peaks <- read.table("./WT_CR_ES_deseq2_sig.bed")
edgeR_peaks <- read.table("./WT_CR_ES_edgeR_sig.bed")

# 给列命名（如果你的 bed 文件确实是这 5 列）
colnames(deseq2_peaks) <- c("chr","start","end","strand","log2FC")

# 2. 转换为 GRanges
gr_deseq2 <- GRanges(seqnames = deseq2_peaks$chr,
                     ranges   = IRanges(start = deseq2_peaks$start,
                                        end   = deseq2_peaks$end),
                     strand   = deseq2_peaks$strand,
                     log2FC   = deseq2_peaks$log2FC)

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


# 筛选距离 TSS 小于 2kb 的 peaks
peakAnno_nearTSS <- subset(peakAnno_df, abs(distanceToTSS) < 3000)

peakAnno_nearTSS_sorted <- peakAnno_nearTSS[order(-abs(peakAnno_nearTSS$log2FC)), ]

# 查看前几行
head(peakAnno_nearTSS_sorted)

write.csv(peakAnno_nearTSS_sorted,
            file = "./WT_CR_ES_deseq2_sig_with_cuttag_nearTSS3kb_sorted.csv",
            quote = FALSE, col.names = TRUE, row.names = FALSE)
