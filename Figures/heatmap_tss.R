# Find peaks in each S phage
gff<-read.table("~/Desktop/Rfiles/idr_peaks/TSS_RT.bed",comment.char = "")
colnames(gff)[c(1,2,3,7)] <- c("chr","start","end","RT")

ES <- gff %>%
  filter(RT == "E") %>%
  select(c(1,2,3,4,5,6))

EMS <- gff %>%
  filter(RT == "EM") %>%
  select(c(1,2,3,4,5,6))

MS <- gff %>%
  filter(RT == "M") %>%
  select(c(1,2,3,4,5,6))

MLS <- gff %>%
  filter(RT == "ML") %>%
  select(c(1,2,3,4,5,6))

LS <- gff %>%
  filter(RT == "L") %>%
  select(c(1,2,3,4,5,6))


write.table(ES,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/ES_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(EMS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/EMS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MLS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MLS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(LS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/LS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)