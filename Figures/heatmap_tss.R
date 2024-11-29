# Find peaks in each S phage
gff<-read.table("~/Desktop/Rfiles/idr_peaks/TSS_RT.bed",comment.char = "")

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