library(tidyverse)
library(ggplot2)

# annotate the genome
positions <- c(0,7,14,19,24,30,44,53,64,92,109,117,124,132,141)
widths <- diff(positions)
transposon <- read.fwf(file="~/Desktop/Rfiles/Rice_MSU7.fasta.std6.9.5.out" , widths=widths)
transposon_rmspace <- as.data.frame(apply(transposon, 2, function(x) gsub(" ","",x)))

transposon_rmspace <- transposon_rmspace %>%
  set_names(transposon_rmspace[2,]) %>% 
  select(c(5,6,7,10)) %>%
  slice(-c(1,2,3)) %>%  
  set_names(c("chr","start","end","family"))

transposon_rmspace[,1] <- gsub("Chr","",transposon_rmspace[,1])
formatC(transposon_rmspace[,1], flag = '0', width = 2)
transposon_rmspace[,1] <- unlist(lapply(as.numeric(transposon_rmspace[,1]),function(x) formatC(x,flag = '0', width = 2)))
transposon_rmspace[,1] <- paste("chr",transposon_rmspace[,1],sep = "")
transposon_rmspace[,5] <- "transposon"
write.table(transposon_rmspace,"~/Desktop/Rfiles/transposon_rmspace.bed",row.names = F,quote=F,sep = "\t",col.names = F)

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
  select(chr,start,end,V7,family) %>%
  set_names(c("chr","start","end","RT","type"))

class1 <- c("LTR/Copia","LTR/Gypsy","LTR/Solo","LTR/TRIM","LTR/unknown",
            "LINE/unknown","SINE/unknown")

transposon_rmspace <- transposon_rmspace %>%
  mutate(type=case_when(type %in% class1 ~ "Class I elements",
                        TRUE ~ "Class II elements"   ))

transposon_class <- transposon_rmspace %>%
  select(chr,start,end,type)
annotate <- rbind(transposon_class,TSS)
# write.table(annotate,"~/Desktop/Rfiles/peak_unit/annotate.bed",row.names = F,quote=F,sep = "\t",col.names = F)

# merge annotate and unannotate file
annotate <- read.table("~/Desktop/Rfiles/peak_unit/annotate_RT.bed", sep = "\t")
annotate <- annotate %>%
  setNames(c("chr","start","end","RT","feature"))

unannotate <- read.table("~/Desktop/Rfiles/peak_unit/unannotated_RT.bed")
unannotate$feature <- "unannotate"
unannotate <- unannotate %>%
  mutate(feature = "unannotate") %>%
  setNames(c("chr","start","end","RT","feature"))

feature_RT <- rbind(unannotate,annotate)
feature_RT <- feature_RT %>%
  arrange(across(everything())) %>%
  mutate(feature =case_when( feature == "gene" ~ "Gene",
                             feature == "Class I elements" ~ "Class I elements",
                             feature == "Class II elements" ~ "Class II elements",
                             feature == "unannotate" ~ "Unannotate")) 

# peaks_reads or peaks_reads_ZH11
peaks_reads <- read.table("~/Desktop/Rfiles/idr_peaks/peaks_reads.txt", header = T)
feature_RT <- feature_RT %>%
  arrange(across(everything()))

number_peaks_ES <- length(which(peaks_reads$RT == "E"))
number_peaks_EMS <- length(which(peaks_reads$RT == "EM"))
number_peaks_MS <- length(which(peaks_reads$RT == "M"))
number_peaks_MLS <- length(which(peaks_reads$RT == "ML"))
number_peaks_LS <- length(which(peaks_reads$RT == "L"))
number_peaks_ELS <- length(which(peaks_reads$RT == "EL"))
number_peaks_EMLS <- length(which(peaks_reads$RT == "EML"))

total_number<- number_peaks_ES + number_peaks_EMS + number_peaks_MS + number_peaks_MLS +
  number_peaks_LS + number_peaks_ELS + number_peaks_EMLS

feature_RT_freq<- as.data.frame(table(feature_RT$feature,feature_RT$RT))
colnames(feature_RT_freq)[1:2] <- c("feature","RT")

RT_freq<- feature_RT_freq %>% group_by(feature)%>%
  summarise(percent = Freq/sum(Freq), RT = RT )

RT_freq$RT = factor(RT_freq$RT,levels = c("E","EM","M","ML","L"))
RT_freq <- RT_freq %>%
  filter(!is.na(RT)) %>%
  mutate(percent.fix = case_when(RT == "E" ~ percent * total_number /number_peaks_ES,
                                 RT == "EM" ~ percent * total_number /number_peaks_EMS,
                                 RT == "M" ~ percent * total_number /number_peaks_MS,
                                 RT == "ML" ~ percent * total_number /number_peaks_MLS,
                                 RT == "L" ~ percent * total_number /number_peaks_LS))

# plot
modify_Figure<- RT_freq %>%
  ggplot(aes(x = feature, y = percent.fix, fill = RT))+
  geom_bar( stat = "identity",colour = "black",position = "dodge")+
  scale_fill_manual(values=c(E = "#2250F1", EM = "#28C5CC", M = "#1A8A12" , ML = "#FFFD33", L = "#FB0018", EL = "#FFEDA0", EML = "#FAB427"))+
  theme_classic() +ggtitle("Replication Times for Chromatin-Related Features")+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Feature")+ylab("Segment coverage")
modify_Figure
# ggsave("~/Desktop/photo//NIP_RT7.pdf", modify_Figure , width = 8, height = 5)      