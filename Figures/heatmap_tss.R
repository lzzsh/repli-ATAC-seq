# Find peaks in each S phage
gff<-read.table("~/Desktop/Rfiles/idr_peaks/TSS_RT.bed",comment.char = "")
colnames(gff)[c(1,2,3,7)] <- c("chr","start","end","RT")

ES <- gff %>%
  filter(RT == "E") 

EMS <- gff %>%
  filter(RT == "EM") %>%

MS <- gff %>%
  filter(RT == "M") %>%

MLS <- gff %>%
  filter(RT == "ML") %>%

LS <- gff %>%
  filter(RT == "L") %>%


write.table(ES,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/ES_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)

write.table(EMS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/EMS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MLS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MLS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(LS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/LS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)