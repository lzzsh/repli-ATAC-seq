library(ChIPseeker)
library(ggplot2)
files <- getSampleFiles()
files <- list.files("~/Desktop/R所需文件/peaks/")
setwd("~/Desktop/R所需文件/peaks")
peak=GenomicRanges::GRangesList(ES=readPeakFile(files[4]),MS=readPeakFile(files[5]),LS=readPeakFile(files[6]))
covplot(peak, weightCol="V5",chrs=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12"))+ 
  facet_grid(chr ~ .id)
library(ggimage)
library(ggupset)
library(GenomicFeatures)
ricegff<-makeTxDbFromGFF("F:/毕设/Arabidopsis/all_DIY.gff3")
peakAnno1 <- annotatePeak(files[[4]], 
                         tssRegion=c(-3000, 3000),
                         TxDb=ricegff, annoDb=ricegff)
plotAnnoPie(peakAnno1)
peakAnno2 <- annotatePeak(files[[5]], 
                         tssRegion=c(-3000, 3000),
                         TxDb=ricegff, annoDb=ricegff)
plotAnnoPie(peakAnno2)
peakAnno3 <- annotatePeak(files[[6]], 
                         tssRegion=c(-3000, 3000),
                         TxDb=ricegff, annoDb=ricegff)
plotAnnoPie(peakAnno3)




library(ggpubr)
annos<-ggarrange(anno3, anno1, anno2, 
          labels = c("ES", "MS", "LS"),
          ncol = 2, nrow = 2)


peakAnnoList <- lapply(files, annotatePeak, TxDb=ricegff,tssRegion=c(-3000, 3000), verbose=FALSE)
for (i in 1:length(peakAnno1)){
  names(data)[i] <- paste("data_", i, sep = "")
}
names(peakAnnoList) <- c("Mid S","G2","G1","Mid S(EdU)","Late S","Early S")
peakAnnoList<-peakAnnoList[c("G1","Mid S","G2","Early S","Mid S(EdU)","Late S")]
png5<-plotDistToTSS(peakAnnoList,title = "Distribution of peaks relative to TSS")
# ggsave("5.png",png5, width = 8, height = 5, dpi = 600)
