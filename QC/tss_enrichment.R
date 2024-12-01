library(GenomicRanges)
library(GenomicFeatures)
setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/tss_plot")
gff3 <- read.table('/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/all_DIY.gff3',stringsAsFactors=F,sep='\t')
ids <- gff3$V3 == 'gene'
genes <- data.frame(chr=gff3$V1[ids],start=gff3$V4[ids],end=gff3$V5[ids],strand=gff3$V7[ids],featrure=gff3$V3[ids],stringsAsFactors=F)
genes.gr <- makeGRangesFromDataFrame(genes,keep.extra.columns=T)
promoters.gr <-promoters(genes.gr,upstream=0,downstream=1)

chr_quant <- function(x){
data_name <- basename(x)
data_name <- gsub('.bam','',data_name)
comm <- sprintf('samtools view -bh -F 260 %s |python ../scripts/extract_quan_in_flag_260_chr_pos.py',x)
data<-matrix(scan(fa<-pipe(comm),''),ncol=2,byrow=T)
colnames(data)<- c('chr','pos')
data <- data.frame(chr=data[,1],start=as.numeric(data[,2]),end=as.numeric(data[,2]),stringsAsFactors=F)
data.gr <- makeGRangesFromDataFrame(data)
dis <-distanceToNearest(data.gr,promoters.gr,ignore.strand=T)
res <- data.frame(promoters.gr[subjectHits(dis)]@ranges@start,data[queryHits(dis),2],promoters.gr[subjectHits(dis)]@strand,stringsAsFactors=F)
ids1 <- res[,3] == "+"
dis1 <- res[ids1,2] - res[ids1,1] 
ids2 <- res[,3] == "-"
dis2 <- res[ids2,1] - res[ids2,2] 
dis <- c(dis1,dis2)
assign(sprintf('%s_tss',data_name),dis)
save(list=sprintf('%s_tss',data_name),file=sprintf('./%s_tss.RData',data_name))
dis_3k<- dis[dis > -3000 & dis < 3000]
x <- table(dis_3k)
x <- cbind(as.numeric(names(x)),x)
sf <- mean(x[x[,1] < -2500| x[,1] >2500,2])
enrichscore <- max(x[,2]/sf)
write.table(enrichscore,file=paste(data_name,'_enrichscore',sep=''),quote = F,row.names = F,col.names =F)
png(file = paste(data_name,'_tss.png',sep=''))
par(mar=c(5,8,5,5),cex=2)
plot(x[,1],x[,2]/sf,type='l',xlab = "Distance to TSS",ylab = "Enrichment")
dev.off()
}

file_list <- function(file){
#file <- read.table(x,stringsAsFactors=F)
for(i in file){
chr_quant(i)
}
}

#file_list(Args[6])
files <- list.files('/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq_control/bwa_all_rawdata','*_q30.bam',full.names=T)
file_list(files)
~