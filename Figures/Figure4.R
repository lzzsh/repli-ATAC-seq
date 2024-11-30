library(dplyr)
tf_location<-read.table("~/Desktop//Rfiles/motif_class.txt")
tf_location<-tf_location[,c(5,6,7,8,4,9,10,11)]
colnames(tf_location)<-c("Chr","Start","End","TF","RT","Id","Strand","PeakID")
tf_location[,5]<-gsub("S","",tf_location[,5])
tf_location_nondup<-distinct(tf_location,Chr,TF,RT,Id,Strand,PeakID, .keep_all= TRUE)
# write.table(table(tf_location_nondup[,4],tf_location_nondup[,5]),"~/Desktop/Rfiles/tf_location_clssfied.bed",row.names = T,quote=F,sep = "\t",col.names = T)

tf_location_nondup<-as.data.frame(read.table("~/Desktop/Rfiles/tf_location_clssfied.bed",header = T))
tf_location_nondup<-tf_location_nondup[,c("E","EM","M","ML","L","EL","EML")]
tf_location_nondup = t(apply(tf_location_nondup,1,function(x){x/sum(x)}))
# tf_location_nondup <- t(apply(tf_location_nondup,1,function(x){log2(x/t(RT_freq$Percent))}))

tf_family<-read.table("~/Desktop/Rfiles/fimo.bed",header = T ,sep = "\t")
tf_family<-tf_family[,c(2,3)]
tf_family<-distinct(tf_family,motif_alt_id ,TF_family,.keep_all= TRUE)

data.1$motif_alt_id<-c(rownames(data.1))
data.1<-merge(data.1,tf_family,by="motif_alt_id")
data.1<-data.1[,c(2,3,4,5,6,7,8,1,9)]
colnames(data.1)[1:7]<-c("E","EM","M","ML","L","EL","EML")
rownames(data.1)<-data.1$motif_alt_id
data.1<-data.1[,-8]

library(pheatmap)
library(vegan)
data.1 <- as.data.frame(decostand(tf_location_nondup,"standardize",MARGIN = 2)) 

colnames(tf_location_nondup) <- c("E","EM","M","ML","L","EL","EML")
pic_heatmap<-pheatmap(data.1,show_rownames = FALSE,show_colnames = TRUE,display_numbers = matrix(ifelse(data.1 > 2, ""," "), nrow(data.1)))

# ggsave("4.pdf",pic_heatmap, width = 8, height = 5)