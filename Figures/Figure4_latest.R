library(dplyr)
library(pheatmap)
library(vegan)
library(readxl)
library(readxl)
library(tidyverse)
library(org.At.tair.db)

tf_location <- read.table("/Users/lzz/Desktop/Rfiles/motif_RT.bed", sep = "\t")
colnames(tf_location)[8] <- "ID"
symbol <- read.table("/Users/lzz/Desktop/Rfiles/symbol.txt")
colnames(symbol) <- c("ID","TF")
tf_location <- merge(tf_location,symbol,by="ID")

# normalization
peaks_reads <- read.table("~/Desktop/Rfiles/idr_peaks/ZH11_RT.gff3", sep = "\t")
colnames(peaks_reads) <- c("chr","start","end","RT")
RT_freq <- peaks_reads %>% group_by(RT)%>%
  summarise(Sum = sum(end) - sum(start))%>%
  filter(!(RT %in% c("EL","EML","ML")))%>%
  mutate(percent = Sum/sum(Sum), RT = factor(RT,levels = c("E", "EM", "M","L")))%>%
  filter(!(is.na(RT)))

tf_location <- tf_location %>%
  mutate(combination = paste(V2, V3, sep = "_")) %>%  
  group_by(combination) %>%                                
  mutate(peak = cur_group_id()) %>%                       
  ungroup() %>%                                          
  mutate(peak = paste0("peak", dense_rank(peak))) %>%      
  dplyr::select(-combination)                                     

tf_location <- tf_location %>%
  dplyr::select(c(6,7,8,14,5,9,1,15)) %>%
  setNames(c("chr","start","end","TF","RT","Strand","ID","PeakID")) 

tf_location_nondup <- distinct(tf_location,chr,TF,PeakID,ID, .keep_all= TRUE) %>%
  distinct(chr,TF,start, .keep_all = TRUE) %>%
  distinct(chr,TF,end, .keep_all = TRUE) %>%
  filter(RT %in% c("E","EM","M","L")) %>%
  arrange(chr, start, end) 

write.table(table(tf_location_nondup$TF,tf_location_nondup$RT),"~/Desktop/Rfiles/peak_unit/tf_location_classfied.bed",row.names = T,quote=F,sep = "\t",col.names = T)

tf_location_nondup <- as.data.frame(read.table("/Users/lzz/Desktop/Rfiles/peak_unit/tf_location_classfied.bed",header = T))
# real_bs <- read_xlsx("~/Desktop/Rfiles/peak_unit/x39_49_wx11_heat.xlsx")
# tf_location_nondup <- rbind(tf_location_nondup,real_bs)
tf_location_nondup <- tf_location_nondup[,c("E","EM","M","L")]
tf_location_nondup <- as.data.frame(t(apply(tf_location_nondup, 1, function(x) x / sum(x))))
# tf_location_nondup = t(apply(tf_location_nondup,1,function(x){x/sum(x)}))
a <- RT_freq$percent
tf_location_nondup[,1] <- tf_location_nondup[,1]/a[1]
tf_location_nondup[,2] <- tf_location_nondup[,2]/a[2]
tf_location_nondup[,3] <- tf_location_nondup[,3]/a[4]
tf_location_nondup[,4] <- tf_location_nondup[,4]/a[3]

# tf_location_nondup <- apply(tf_location_nondup, 2, function(x, percent) x / percent, percent = c(a[1],a[2],a[4],a[5],a[3]))

tf_family <- read.table("/Users/lzz/Desktop/Rfiles/fimo.bed",header = T ,sep = "\t")
tf_family <- tf_family[,c(2,3)]
tf_family <- distinct(tf_family,motif_alt_id ,TF_family,.keep_all= TRUE)


# plot
data.1 <- as.data.frame(decostand(tf_location_nondup,"standardize",MARGIN = 2)) 
pic_heatmap<-pheatmap(data.1,show_rownames = FALSE,show_colnames = TRUE,
                      display_numbers = matrix(ifelse(data.1 > 2, ""," "), 
                                               nrow(data.1)),cluster_rows = FALSE)
#,color = colorRampPalette(c("#547297", "#8C9EBA", "#D9E0E7","#F3DBD6","#DA8F87"
#                           ,"#D54846"))(100))
ggsave("~/Documents/GitHub/repli-ATAC-seq/output/Figures/tf_location_ZH11_4RT.pdf", pic_heatmap , width = 8, height = 5)

rownames(data.1)[rownames(data.1) == "ATHB-40"] <- "HB-5"
data.1$motif_alt_id<-c(rownames(data.1))
data.1<-merge(data.1,tf_family,by="motif_alt_id", all.x=TRUE)
data.1<-data.1[,c(2,3,4,5,1,6)]
colnames(data.1)[1:4]<-c("E","EM","M","L")
rownames(data.1)<-data.1$motif_alt_id
# write.table(data.1,"~/Desktop/Rfiles/peak_unit/tf_location_NIP_4RT.csv",quote=F,sep = ",")

# annotation of tf_location
trans<-read.csv("~/Desktop/repli-ATAC-seq/rice2ara.csv")
geneid<-select(org.At.tair.db, keys = c(t(trans[,2])) ,
               column = c('ENTREZID', 'SYMBOL', 'REFSEQ'), keytype = 'TAIR')
geneid<-distinct(geneid,TAIR,SYMBOL)
colnames(data.1)[5] <- "SYMBOL"
a<-merge(data.1, geneid ,by = "SYMBOL" , all.x = TRUE)
matched <- str_subset(data.1$SYMBOL, "^AT[0-9]G[0-9]{5}$")
a[which(a$SYMBOL %in% matched),7] <- matched
colnames(trans)[2]<-"TAIR"
a<-merge(a,trans,by = "TAIR", all.x = TRUE)
msu_to_symbol<-read_xlsx("~/Desktop/Rfiles/gene_function.xlsx")
msu_to_symbol<-msu_to_symbol[,c(1,8)]
colnames(a)[8]<-"MSU"
a<-merge(a,msu_to_symbol,by = "MSU", all.x = TRUE)
a <- a %>%
  dplyr::select(c(1,2,3,8,4,5,6,7,9)) %>%
  setNames(c("MSU","TAIR","SYMBOL_TAIR","TF_family_TAIR","E","EM","M","L","SYMBOL_NIP"))

# add gene expression information
rpkm <- read.table("~/Desktop/Rfiles/peak_unit/NIP_rep1_rep2_FPKM.featureCounts.matrix", skip = 186)
rpkm <- rpkm[,c(1,3,5)]
colnames(rpkm) <- c("MSU","FPKM.1","FPKM.2") 
a <- merge(a,rpkm,by="MSU",all.x=TRUE)
# write.table(a,"~/Documents/GitHub/repli-ATAC-seq/output/tf_location_NIP_4RT_all.csv",quote=F,sep = ",",row.names = F)
