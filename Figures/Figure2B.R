library(dplyr)
library(ggplot2)

# load histone_RT.bed file 
peaks_reads <- read.table("~/Desktop/Rfiles/idr_peaks/ZH11_RT.gff3", sep = "\t")
histone <- read.table("~/Desktop/Rfiles/peak_unit/Histone_RT.bed", sep = "\t")
colnames(peaks_reads) <- c("chr","start","end","RT")
colnames(histone)[1:3]<- c("chr","start","end")

histone <- merge(histone,peaks_reads[,c("chr","start","end","RT")],by=c("chr","start","end")) 
histone <- histone %>%
  dplyr::select(5,6,7,8,10)%>%
  filter(!(RT %in% c("EL","EML")))%>%
  setNames(c("chr","start","end","Modify","RT")) %>%
  arrange(across(everything()))

# normalization
RT_freq <- peaks_reads %>% group_by(RT)%>%
  summarise(Sum = sum(end) - sum(start))%>%
  filter(!(RT %in% c("EL","EML")))%>%
  mutate(percent = Sum/sum(Sum), RT = factor(RT,levels = c("E", "EM", "M","ML","L")))

histone_freq <- as.data.frame(table(histone$Modify,histone$RT)) %>%
  setNames(c("Modify","RT","Freq")) %>%
  group_by(Modify) %>%
  reframe(RT = RT, percent = Freq / sum(Freq))

modify_plot <- merge(histone_freq, RT_freq, by = "RT") %>%
  dplyr::select(1,2,3,5) %>%
  mutate(percent.fix = percent.x / percent.y)

modify_plot$RT = factor(modify_plot$RT,levels = c("E", "EM", "M","ML","L"))
modify_plot$Modify = factor(modify_plot$Modify,levels = c("H3K4me1","H3K27ac","H3K4me3",
                                                          "H3K27me3","H3K9me2"))

# plot
modify_Figure<- modify_plot %>%
  ggplot(aes(x = Modify, y = percent.fix, fill = RT))+
  geom_bar( stat = "identity",colour = "black",position = "dodge")+
  scale_fill_manual(values=c(E = "#2250F1", EM = "#28C5CC", M = "#1A8A12" , ML = "#FFFD33", L = "#FB0018", EL = "#FFEDA0", EML = "#FAB427"))+
  theme_classic() +ggtitle("Replication Times for Chromatin-Related Features")+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Histone mark signature")+ylab("% of RT class with mark signature")
modify_Figure                                  
ggsave("~/Documents/Github/repli-ATAC-seq/output/Figures/Figure2B.pdf", modify_Figure , width = 8, height = 5)
