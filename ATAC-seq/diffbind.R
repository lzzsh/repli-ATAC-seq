setwd("/storage2/liuxiaodongLab/liaozizhuo/Projects/ATAC-seq-CR/diffbind/")
library(DiffBind)
# 1. 读入 sampleSheet
samplesheet <- read.csv("./sampleSheet1.csv")
samplesheet <- samplesheet[c(4,8),]

dbObj <- dba(sampleSheet = samplesheet)

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dba.plotPCA(dbObj, attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)
dba.plotVenn(dbObj, dbObj$masks$All)

# Establishing a contrast
dbObj <- dba.contrast(dbObj, categories=c(DBA_CONDITION), minMembers=1)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
# summary of results
dba.show(dbObj, bContrasts=T)
# overlapping peaks identified by the two different tools (DESeq2 and edgeR)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

## save results
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="./WT_CR_MS_edgeR_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="./WT_CR_MS_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)



# 合并成 matrix
mat <- dba.peakset(dbObj, bRetrieve=TRUE)

# 转成 data.frame
mat_df <- as.data.frame(mat)

# 比如你自己加了 log2FC 列
mat_df$log2FC <- log2((mat$WT.LS.1 + 1)/(mat$TCX2_1.LS.1 + 1))

# 用 dplyr 过滤
mat_df_filtered <- mat_df %>%
  filter(abs(log2FC) > 1)

