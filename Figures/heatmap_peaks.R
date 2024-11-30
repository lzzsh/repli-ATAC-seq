# Find peaks in each S phage
gff<-read.table("~/Desktop/Rfiles/idr_peaks/ZH11_gff.gff3",skip = 2,comment.char = "")
gff<-gff[,c(1,4,5,9)]
colnames(gff)<-c("chr","start","end","RT")
gff$RT<-lapply(gff$RT, function(x) unlist(strsplit(x,";"))[1])
gff$RT<-lapply(gff$RT, function(x) unlist(strsplit(x,"="))[2])
gff<-gff[ !(gff$chr %in% c("chrUn","chrSy")),] # remove chrSy and chrUn
gff[,4]<-gsub("S","",gff[,4])

ES <- gff %>%
  filter(RT == "E") %>%
  select(chr,start,end)
EMS <- gff %>%
  filter(RT == "EM") %>%
  select(chr,start,end)
MS <- gff %>%
  filter(RT == "M") %>%
  select(chr,start,end)
MLS <- gff %>%
  filter(RT == "ML") %>%
  select(chr,start,end)
LS <- gff %>%
  filter(RT == "L") %>%
  select(chr,start,end)

write.table(ES,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/ES_peaks.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(EMS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/EMS_peaks.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MS_peaks.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MLS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MLS_peaks.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(LS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/LS_peaks.txt",row.names = F,quote=F,sep = "\t",col.names = F)