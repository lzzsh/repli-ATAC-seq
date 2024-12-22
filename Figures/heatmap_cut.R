# x39
library(dplyr)
gff<-read.table("../x39_cut_out/x39_RT_top3000.bed",comment.char = "")
colnames(gff)<-c("chr","start","end","RT")

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

write.table(ES,"../x39_cut_out/x39_ES_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(EMS,"../x39_cut_out/x39_EMS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MS,"../x39_cut_out/x39_MS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MLS,"../x39_cut_out/x39_MLS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(LS,"../x39_cut_out/x39_LS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)

# x49
gff<-read.table("../x49_cut_out/x49_RT_top3000.bed",comment.char = "")
colnames(gff)<-c("chr","start","end","RT")

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

write.table(ES,"../x49_cut_out/x49_ES_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(EMS,"../x49_cut_out/x49_EMS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MS,"../x49_cut_out/x49_MS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MLS,"../x49_cut_out/x49_MLS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(LS,"../x49_cut_out/x49_LS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)

# xw11
gff<-read.table("../xw11_cut_out/xw11_RT_top3000.bed",comment.char = "")
colnames(gff)<-c("chr","start","end","RT")

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

write.table(ES,"../xw11_cut_out/xw11_ES_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(EMS,"../xw11_cut_out/xw11_EMS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MS,"../xw11_cut_out/xw11_MS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(MLS,"../xw11_cut_out/xw11_MLS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)
write.table(LS,"../xw11_cut_out/xw11_LS_peaks_top3000.txt",row.names = F,quote=F,sep = "\t",col.names = F)