# Find peaks in each S phage

# IDR_peaks
library(dplyr)
gff<-read.table("~/Desktop/Rfiles/idr_peaks/TSS_RT.bed",comment.char = "")
colnames(gff)[c(1,2,3,7)] <- c("chr","start","end","RT")

ES <- gff %>%
  filter(RT == "E") %>%
  dplyr::select(c(1,2,3,4,5,6))

EMS <- gff %>%
  filter(RT == "EM") %>%
  dplyr::select(c(1,2,3,4,5,6))

MS <- gff %>%
  filter(RT == "M") %>%
  dplyr::select(c(1,2,3,4,5,6))

MLS <- gff %>%
  filter(RT == "ML") %>%
  dplyr::select(c(1,2,3,4,5,6))

LS <- gff %>%
  filter(RT == "L") %>%
  dplyr::select(c(1,2,3,4,5,6))


write.table(ES,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/ES_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(EMS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/EMS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MLS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MLS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(LS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/LS_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)

# Control_IDR_peaks
library(dplyr)
gff<-read.table("~/Desktop/Rfiles/idr_peaks/TSS_RT_control.bed",comment.char = "")
colnames(gff)[c(1,2,3,7)] <- c("chr","start","end","RT")

ES <- gff %>%
  filter(RT == "E") %>%
  dplyr::select(c(1,2,3,4,5,6))

EMS <- gff %>%
  filter(RT == "EM") %>%
  dplyr::select(c(1,2,3,4,5,6))

MS <- gff %>%
  filter(RT == "M") %>%
  dplyr::select(c(1,2,3,4,5,6))

MLS <- gff %>%
  filter(RT == "ML") %>%
  dplyr::select(c(1,2,3,4,5,6))

LS <- gff %>%
  filter(RT == "L") %>%
  dplyr::select(c(1,2,3,4,5,6))


write.table(ES,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/ES_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(EMS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/EMS_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MS_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MLS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MLS_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(LS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/LS_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)

# ZH11-2_peaks
library(dplyr)
gff<-read.table("~/Desktop/Rfiles/idr_peaks/TSS_RT_2_control.bed",comment.char = "")
colnames(gff)[c(1,2,3,7)] <- c("chr","start","end","RT")

ES <- gff %>%
  filter(RT == "E") %>%
  dplyr::select(c("chr","start","end","V4","V5","V6"))

EMS <- gff %>%
  filter(RT == "EM") %>%
  dplyr::select(c(1,2,3,4,5,6))

MS <- gff %>%
  filter(RT == "M") %>%
  dplyr::select(c(1,2,3,4,5,6))

MLS <- gff %>%
  filter(RT == "ML") %>%
  dplyr::select(c(1,2,3,4,5,6))

LS <- gff %>%
  filter(RT == "L") %>%
  dplyr::select(c(1,2,3,4,5,6))


write.table(ES,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/ES_2_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(EMS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/EMS_2_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MS_2_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MLS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/MLS_2_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(LS,"~/Documents/GitHub/repli-ATAC-seq/output/heatmap/LS_2_control_tss.txt",row.names = F,quote=F,sep = "\t",col.names = F)