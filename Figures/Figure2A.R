library(tidyverse)
annotate <- read.table("~/Desktop/R所需文件/peaks划分/annotate_RT.bed")
annotate <- annotate %>%
  select(c(1,2,3,7)) %>%
  setNames(c("chr","start","end","feature"))

unannotate <- read.table("~/Desktop/R所需文件/peaks划分/unannotate.bed")
unannotate <- unannotate %>%
  mutate(feature = "unannotate") %>%
  setNames(c("chr","start","end","feature"))

feature <- rbind(unannotate,annotate)
feature <- feature %>%
  arrange(across(everything())) %>%
  mutate(feature =case_when( feature == "gene" ~ "Gene",
                             feature == "class1" ~ "Class I elements",
                             feature == "class2" ~ "Class II elements",
                             feature == "unannotate" ~ "Unannotate")) 



#peaks_reads or peaks_reads_ZH11
feature_RT<-merge(feature,peaks_reads[,c("chr","start","end","RT")],by=c("chr","start","end")) 
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
  mutate(percent.fix = case_when(RT == "E" ~ percent * total_number /number_peaks_ES,
                                 RT == "EM" ~ percent * total_number /number_peaks_EMS,
                                 RT == "M" ~ percent * total_number /number_peaks_MS,
                                 RT == "ML" ~ percent * total_number /number_peaks_MLS,
                                 RT == "L" ~ percent * total_number /number_peaks_LS))

modify_Figure<- RT_freq %>%
  ggplot(aes(x = feature, y = percent.fix, fill = RT))+
  geom_bar( stat = "identity",colour = "black",position = "dodge")+
  scale_fill_manual(values=c(E = "#2250F1", EM = "#28C5CC", M = "#1A8A12" , ML = "#FFFD33", L = "#FB0018", EL = "#FFEDA0", EML = "#FAB427"))+
  theme_classic() +ggtitle("Replication Times for Chromatin-Related Features")+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Feature")+ylab("Segment coverage")
modify_Figure