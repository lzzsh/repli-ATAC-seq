library(ChIPseeker)
library(ggplot2)
library(ggimage)
library(ggupset)
library(GenomicFeatures)
library(ggpubr)

options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = F)

# load 6 files
files <- list.files("~/Desktop/Rfiles/peaks/")
setwd("~/Desktop/Rfiles/peaks")

ricegff<-makeTxDbFromGFF("~/Desktop/Rfiles/all_DIY.gff3")

# NIP
peakAnno1 <- annotatePeak(files[[1]], 
                          tssRegion=c(-3000, 3000),
                          TxDb=ricegff, annoDb=ricegff,
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
plotAnnoPie(peakAnno1)
peakAnno2 <- annotatePeak(files[[3]], 
                          tssRegion=c(-3000, 3000),
                          TxDb=ricegff, annoDb=ricegff,
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
plotAnnoPie(peakAnno2)
peakAnno3 <- annotatePeak(files[[2]], 
                          tssRegion=c(-3000, 3000),
                          TxDb=ricegff, annoDb=ricegff,
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
plotAnnoPie(peakAnno3)

# ZH11
peakAnno1 <- annotatePeak(files[[4]], 
                          tssRegion=c(-3000, 3000),
                          TxDb=ricegff, annoDb=ricegff,
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
plotAnnoPie(peakAnno1)

peakAnno2 <- annotatePeak(files[[6]], 
                          tssRegion=c(-3000, 3000),
                          TxDb=ricegff, annoDb=ricegff,
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
plotAnnoPie(peakAnno2)
peakAnno3 <- annotatePeak(files[[5]], 
                          tssRegion=c(-3000, 3000),
                          TxDb=ricegff, annoDb=ricegff,
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
plotAnnoPie(peakAnno3)

# bar plot
peakAnnoList <- lapply(files, annotatePeak, TxDb=ricegff,tssRegion=c(-3000, 3000), verbose=FALSE)
names(peakAnnoList) <- c("NIP_ES","NIP_MS","NIP_LS","Early S","Late S","Mid S(EdU)","G1","G2","Mid S")
peakAnnoList<-peakAnnoList[c("G1","Mid S","G2","Early S","Mid S(EdU)","Late S")]
png5<-plotDistToTSS(peakAnnoList,title = "Distribution of peaks relative to TSS")
png5
# ggsave("5.png",png5, width = 8, height = 5, dpi = 600)