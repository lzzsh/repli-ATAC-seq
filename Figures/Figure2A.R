library(tidyverse)
library(ggplot2)

# annotate the genome
positions <- c(0,7,14,19,24,30,44,53,64,92,109,117,124,132,141)
widths <- diff(positions)
transposon <- read.fwf(file="~/Desktop/Rfiles/Rice_MSU7.fasta.std6.9.5.out" , widths=widths)
transposon_rmspace <- as.data.frame(apply(transposon, 2, function(x) gsub(" ","",x)))

transposon_rmspace <- transposon_rmspace %>%
  set_names(transposon_rmspace[2,]) %>% 
  dplyr::select(c(5,6,7,10)) %>%
  dplyr::slice(-c(1,2,3)) %>%  
  set_names(c("chr","start","end","family"))

transposon_rmspace[,1] <- gsub("Chr","",transposon_rmspace[,1])
formatC(transposon_rmspace[,1], flag = '0', width = 2)
transposon_rmspace[,1] <- unlist(lapply(as.numeric(transposon_rmspace[,1]),function(x) formatC(x,flag = '0', width = 2)))
transposon_rmspace[,1] <- paste("chr",transposon_rmspace[,1],sep = "")
transposon_rmspace[,5] <- "transposon"
# write.table(transposon_rmspace,"~/Desktop/Rfiles/transposon_rmspace.bed",row.names = F,quote=F,sep = "\t",col.names = F)

# merge tss location and transposon
TSS <- read.table("~/Desktop/Rfiles/peak_unit/TSS_location.bed")
TSS[,4] <- "gene"
colnames(TSS) <- c("chr","start","end","type")

transposon_rmspace <- read.table("~/Desktop/Rfiles/transposon_rmspace.bed")
colnames(transposon_rmspace) <- c("chr","start","end","type","family")

a <- read.table("~/Desktop/Rfiles/transposon_RT.bed",sep = "\t")
colnames(a)[1:3] <- c("chr","start","end")

transposon_RT <- merge(a,transposon_rmspace,by=c("chr","start","end")) %>%
  arrange(chr,start,end) %>%
  dplyr::select(chr,start,end,V7,family) %>%
  set_names(c("chr","start","end","RT","type"))

class1 <- c("LTR/Copia","LTR/Gypsy","LTR/Solo","LTR/TRIM","LTR/unknown",
            "LINE/unknown","SINE/unknown")

transposon_rmspace <- transposon_rmspace %>%
  mutate(type=case_when(type %in% class1 ~ "Class I elements",
                        TRUE ~ "Class II elements"   ))

transposon_class <- transposon_rmspace %>%
  dplyr::select(chr,start,end,type)
annotate <- rbind(transposon_class,TSS)
# write.table(annotate,"~/Desktop/Rfiles/peak_unit/annotate.bed",row.names = F,quote=F,sep = "\t",col.names = F)

# merge annotate and unannotate file
annotate <- read.table("~/Desktop/Rfiles/peak_unit/annotate_RT_org.bed", sep = "\t")
annotate <- annotate %>%
  setNames(c("chr","start","end","RT","feature"))

unannotate <- read.table("~/Desktop/Rfiles/peak_unit/unannotated_RT_org.bed")
unannotate$feature <- "unannotate"
unannotate <- unannotate %>%
  mutate(feature = "unannotate") %>%
  setNames(c("chr","start","end","RT","feature"))

feature_RT <- rbind(unannotate,annotate)
feature_RT <- feature_RT %>%
  filter(RT %in% c("ES","ESMS","MS","MSLS","LS")) %>%
  arrange(across(everything())) %>%
  mutate(feature =case_when( feature == "gene" ~ "Gene",
                             feature == "Class I elements" ~ "Class I elements",
                             feature == "Class II elements" ~ "Class II elements",
                             feature == "unannotate" ~ "Unannotate")) 

# peaks_reads or peaks_reads_ZH11
peaks_reads <- read.table("~/Desktop/Rfiles/peak_unit/ZH11_RT_org.gff3", sep = "\t")
colnames(peaks_reads) <- c("chr","start","end","RT")
feature_RT <- feature_RT %>%
  arrange(across(everything()))

# Calculate RT frequency based on genomic length (same as Figure 2B and 2CD)
RT_freq_genome <- peaks_reads %>% group_by(RT) %>%
  summarise(Sum = sum(end) - sum(start)) %>%
  filter(!(RT %in% c("ESLS","ESMSLS"))) %>%
  mutate(percent = Sum/sum(Sum), 
         RT = factor(RT, levels = c("ES", "ESMS", "MS","MSLS","LS")))

feature_RT_freq <- as.data.frame(table(feature_RT$feature,feature_RT$RT))
colnames(feature_RT_freq)[1:2] <- c("feature","RT")

# Calculate feature frequency
feature_freq <- feature_RT_freq %>% group_by(feature) %>%
  summarise(percent = Freq / sum(Freq), RT = RT)

# Normalization (same method as Figure 2B and 2CD)
RT_freq <- merge(feature_freq, RT_freq_genome, by = "RT") %>%
  dplyr::select(1,2,3,5) %>%
  mutate(percent.fix = percent.x / percent.y)  # Normalize by RT frequency

RT_freq$RT = factor(RT_freq$RT,levels = c("ES","ESMS","MS","MSLS","LS"))
RT_freq <- RT_freq %>%
  filter(!is.na(RT))

# plot
Figure2A <- RT_freq %>%
  ggplot(aes(x = feature, y = percent.fix, fill = RT))+
  geom_bar( stat = "identity",colour = "black",position = "dodge")+
  scale_fill_manual(values=c(ES = "#2C5F9E", ESMS = "#68A0D8", MS = "#95BE6C", MSLS = "#E4B660", LS = "#E68364", ESLS = "#B784A7", ESMSLS = "#9B7EB3"))+
  theme_classic() +ggtitle("Replication Times for Chromatin-Related Features")+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Feature")+ylab("% of RT class with feature")
Figure2A
# ggsave("~/Documents/Github/repli-ATAC-seq/output/Figures/Figure2A.pdf", Figure2A , width = 8, height = 5)
