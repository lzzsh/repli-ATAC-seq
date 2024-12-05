library(rtracklayer)
library(dplyr)
library(ggplot2)

gtf_data = import('~/Desktop/Rfiles/all_DIY.gff3')
gtf_data = as.data.frame(gtf_data)
gene_data <- gtf_data[which(!is.na(gtf_data[,11]) & gtf_data[,7]=="gene"),c(1,2,3,5,11)]
colnames(gene_data) <- c("chrom","start","end","strand","Geneid")
gene_data$end <- gene_data$start;gene_data$start <- gene_data$end - 1
# write.table(gene_data,"~/Desktop/Rfiles/gene_data.txt",row.names = F,quote=F,sep = "\t",col.names = F)

TSS_RT <- read.table("~/Desktop/Rfiles/idr_peaks/TSS_RT_2_control.bed") %>%
  dplyr::select(c(1,2,3,6,7)) %>%
  setNames(c("chrom","start","end","strand","RT"))

gene_rpkm <- merge(TSS_RT,gene_data,by=c("chrom","start","end","strand"), all.x=TRUE)
gene_rpkm <- gene_rpkm %>%
  arrange(across(everything()))

rpkm <- read.table("~/Desktop/Rfiles/peak_unit/NIP_rep1_rep2_FPKM.featureCounts.matrix", header = T, skip = 185)
colnames(rpkm) <- c("Geneid","NIP_rep1","FPKM","NIP_rep2","FPKM.1")
rpkm <- rpkm[,c(1,3,5)]
rpkm$FPKM_average <- (rpkm$FPKM + rpkm$FPKM.1) / 2

gene_rpkm<-merge(gene_rpkm[,-8],rpkm,by="Geneid")
gene_rpkm[,10]=0
for ( i in 1:nrow(gene_rpkm) )
{
  if(gene_rpkm[i,9]==0)
  {gene_rpkm[i,10]="F=0"}
  if(gene_rpkm[i,9]>0 & gene_rpkm[i,9]<=1)
  {gene_rpkm[i,10]="0<F<=1"}
  if(gene_rpkm[i,9]>1 & gene_rpkm[i,9]<=10)
  {gene_rpkm[i,10]="1<F<=10"}
  if(gene_rpkm[i,9]>10 & gene_rpkm[i,9]<=100)
  {gene_rpkm[i,10]="10<F<=100"}
  if(gene_rpkm[i,9]>100)
  {gene_rpkm[i,10]="F>100"}
  
}

gene_rpkm<-gene_rpkm[which(gene_rpkm$RT %in% c("E","EM","M","ML","L")),]

rpkm_table<-as.data.frame(table(gene_rpkm[,10],gene_rpkm[,6]))
colnames(rpkm_table)<-c("RPKM","RT","FREQ")
rpkm_table1<-rpkm_table %>% group_by(RPKM)%>%
  summarise(Sum = sum(FREQ))
RPKM_plot<-merge(rpkm_table,rpkm_table1,by="RPKM")
RPKM_plot$Ratio<-RPKM_plot[,3]/RPKM_plot[,4]*100
RPKM_plot<-RPKM_plot[which(RPKM_plot[,2]!="EL" & RPKM_plot[,2]!="EML"),]
RPKM_plot$RT = factor(RPKM_plot$RT,levels = c("E", "EM", "M","ML","L"))
RPKM_plot$RPKM = factor(RPKM_plot$RPKM,levels = c("F=0","0<F<=1","1<F<=10",
                                                  "10<F<=100","F>100"))

# normalization
peaks_reads <- read.table("~/Desktop/Rfiles/idr_peaks/peaks_reads_2_control.txt", header = T, sep = "\t")

RT_freq <- gene_rpkm %>% group_by(RT)%>%
  summarise(Sum = sum(end) - sum(start))%>%
  filter(!(RT %in% c("EL","EML")))%>%
  mutate(percent = Sum/sum(Sum), RT = factor(RT,levels = c("E", "EM", "M","ML","L")))

modify_plot <- merge(RPKM_plot,RT_freq,by="RT")
modify_plot$Ratio <- modify_plot$Ratio/modify_plot$percent/100

# FPKM classfication chart
Figure2C <- modify_plot %>%
  ggplot(aes(x = RPKM, y = Ratio, fill = RT))+
  geom_bar( stat = "identity",colour = "black",position = "dodge")+
  scale_fill_manual(values=c(E = "#2250F1", EM = "#28C5CC", M = "#1A8A12" , ML = "#FFFD33", L = "#FB0018", EL = "#FFEDA0", EML = "#FAB427"))+
  theme_classic() +
  xlab("Gene expression level (FPKM)")+ylab("of genes in expression level")
Figure2C
ggsave("~/Documents/Github/repli-ATAC-seq/output/Figures/Figure2C.pdf", Figure2C , width = 8, height = 5)

# FPKM boxplot
gene_nonzero <- gene_rpkm %>%
  filter(FPKM_average > 0)
gene_nonzero$RT = factor(gene_nonzero$RT,levels = c("E", "EM", "M","ML","L"))
boxplot(FPKM_average ~ RT, data = gene_nonzero, col = c(E = "#2250F1", EM = "#28C5CC", M = "#1A8A12" , ML = "#FFFD33", L = "#FB0018"),outline=FALSE)