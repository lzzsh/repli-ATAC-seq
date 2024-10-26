library(dplyr)
library(readxl)
library(pheatmap)
library(vegan)
library(readxl)
library(tidyverse)
library(org.At.tair.db)

tf_location<-read.table("/Users/lzz/Desktop/R所需文件/peaks划分/tf_motif.bed")
colnames(tf_location)[1:3] <- c("chr","start","end")
tf_location<-merge(tf_location,peaks_reads[,c("chr","start","end","RT")],by=c("chr","start","end")) 

tf_location <- tf_location %>%
  dplyr::select(c(4,5,6,7,8,9,10,11)) %>%
  setNames(c("chr","start","end","TF","Id","Strand","PeakID","RT")) 
  
tf_location_nondup <- distinct(tf_location,chr,TF,Id,PeakID, .keep_all= TRUE) %>%
  distinct(chr,TF,start, .keep_all = TRUE) %>%
  distinct(chr,TF,end, .keep_all = TRUE) %>%
  arrange(chr, start, end)

# write.table(table(tf_location_nondup$TF,tf_location_nondup$RT),"~/Desktop/R所需文件/peaks划分/tf_location_clssfied1.bed",row.names = T,quote=F,sep = "\t",col.names = T)

tf_location_nondup<-as.data.frame(read.table("/Users/lzz/Desktop/R所需文件/peaks划分/tf_location_clssfied1.bed",header = T))
real_bs<-read_xlsx("~/Desktop/R所需文件/peaks划分/x39_49_wx11_heat.xlsx")
tf_location_nondup<-rbind(tf_location_nondup,real_bs)
tf_location_nondup<-tf_location_nondup[,c("E","EM","M","ML","L","EL","EML")]
tf_location_nondup = t(apply(tf_location_nondup,1,function(x){x/sum(x)}))
# tf_location_nondup <- t(apply(tf_location_nondup,2,function(x){log2(x/t(RT_freq$Percent))}))

tf_family<-read.table("/Users/lzz/Desktop/R所需文件/fimo.bed",header = T ,sep = "\t")
tf_family<-tf_family[,c(2,3)]
tf_family<-distinct(tf_family,motif_alt_id ,TF_family,.keep_all= TRUE)




data.1 <- as.data.frame(decostand(tf_location_nondup,"standardize",MARGIN = 2)) 
pic_heatmap<-pheatmap(data.1,show_rownames = FALSE,show_colnames = TRUE,
                      display_numbers = matrix(ifelse(data.1 > 2, ""," "), 
                      nrow(data.1)),cluster_rows = FALSE)
                      #,color = colorRampPalette(c("#547297", "#8C9EBA", "#D9E0E7","#F3DBD6","#DA8F87"
                      #                           ,"#D54846"))(100))
ggsave("~/Desktop/photo/tf_location_ZH11_7RT.png", pic_heatmap , width = 8, height = 5, dpi = 300)

data.1$motif_alt_id<-c(rownames(data.1))
data.1<-merge(data.1,tf_family,by="motif_alt_id")
data.1<-data.1[,c(2,3,4,5,6,7,8,1,9)]
colnames(data.1)[1:7]<-c("E","EM","M","ML","L","EL","EML")
rownames(data.1)<-data.1$motif_alt_id
#write.table(data.1,"~/Desktop/R所需文件/peaks划分/tf_location_NIP_7RT.csv",quote=F,sep = ",")



###########annotation of tf_location#######
trans<-read.csv("~/Desktop/毕设文献/rice2ara.csv")
geneid<-select(org.At.tair.db, keys = c(t(trans[,2])) ,
               column = c('ENTREZID', 'SYMBOL', 'REFSEQ'), keytype = 'TAIR')
geneid<-distinct(geneid,TAIR,SYMBOL)
colnames(data.1)[8] <- "SYMBOL"
a<-merge(data.1, geneid ,by = "SYMBOL" , all.x = TRUE)
matched <- str_subset(data.1$SYMBOL, "^AT[0-9]G[0-9]{5}$")
a[which(a$SYMBOL %in% matched),10] <- matched
colnames(trans)[2]<-"TAIR"
a<-merge(a,trans,by = "TAIR", all.x = TRUE)
msu_to_symbol<-read_xlsx("~/Desktop/R所需文件/gene_function.xlsx")
msu_to_symbol<-msu_to_symbol[,c(1,8)]
colnames(a)[11]<-"MSU"
a<-merge(a,msu_to_symbol,by = "MSU", all.x = TRUE)
a <- a %>%
  dplyr::select(c(1,2,3,11,4,5,6,7,8,9,10,12)) %>%
  setNames(c("MSU","TAIR","SYMBOL_TAIR","TF_family_TAIR","E","EM","M","ML","L","EL","EML","SYMBOL_NIP"))
#write.table(a,"~/Desktop/R所需文件/peaks划分/tf_location_NIP_7RT_all.csv",quote=F,sep = ",",row.names = F)