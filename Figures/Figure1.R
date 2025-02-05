library(dplyr)
library(ggplot2)
library(ggstatsplot)

# load gff which contains the RT location of each open chromatin region
gff <- read.table("~/Desktop/repli-ATAC-seq/segamentation/ZH11.gff3",comment.char = "")
gff <- gff[,c(1,4,5,9)]
colnames(gff)<-c("chr","start","end","RT")
gff$RT <- lapply(gff$RT, function(x) unlist(strsplit(x,";"))[1])
gff$RT <- lapply(gff$RT, function(x) unlist(strsplit(x,"="))[2])
gff <- gff[ !(gff$chr %in% c("chrUn","chrSy")),] 
gff$length <- gff$end-gff$start+1
gff<-gff[!(is.na(gff$RT)),]
gff[,4] <- gsub("S","",gff[,4])

RT_length <- gff %>% group_by(RT) %>%
  summarise(Sum = sum(length))
RT_length$freq <- RT_length[,2]/sum(RT_length[,2])

# pie chart
gff$RT <- factor(gff$RT,levels = c("EML","EL","L","ML","M","EM","E"))
library(ggrepel)

Figure1A <- ggplot(gff, aes(x = "", fill = RT)) +
  geom_bar(width = 1, aes(weight = length), stat = "count") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(E = "#2250F1", EM = "#28C5CC", M = "#1A8A12", ML = "#FFFD33", L = "#FB0018", EL = "#FFEDA0", EML = "#FAB427")) +
  labs(title = "Partition of total genome replication into seven RT segment classes") +
  theme_void()
Figure1A
# ggsave("~/Documents/Github/repli-ATAC-seq/output/Figures/Figure1A.pdf", Figure1A , width = 8, height = 5)

gff$length <- log2((gff$end-gff$start+1))

# boxplot
gff$RT <- factor(gff$RT,levels = c("E","EM","M","ML","L","EL","EML"))
Figure1B <- ggplot(data = gff, aes(x = RT, y = length))+   
  geom_boxplot(aes(color = RT), size = 1.1)+ 
  scale_color_manual(values=c(E = "#2250F1", EM = "#28C5CC", M = "#1A8A12" , ML = "#FFFD33", L = "#FB0018", EL = "#FFEDA0", EML = "#FAB427"))+
  scale_fill_manual(values=c(E = "#2250F1", EM = "#28C5CC", M = "#1A8A12" , ML = "#FFFD33", L = "#FB0018", EL = "#FFEDA0", EML = "#FAB427"))+
  theme_classic()+ ggtitle("Size distribution of replication timing segments.")+theme(plot.title = element_text(hjust = 0.5))+ylab("log2(Length)")
Figure1B   
# ggsave("~/Documents/Github/repli-ATAC-seq/output/Figures/Figure1B.pdf", Figure1B , width = 8, height = 5)
