library(dplyr)
library(pheatmap)
library(vegan)
library(readxl)
library(readxl)
library(tidyverse)
library(org.At.tair.db)

tf_location <- read.table("/Users/lzz/Desktop/Rfiles/motif_RT_org.bed", sep = "\t")
colnames(tf_location)[8] <- "ID"
symbol <- read.table("/Users/lzz/Desktop/Rfiles/symbol.txt")
colnames(symbol) <- c("ID","TF")
tf_location <- merge(tf_location,symbol,by="ID")

# normalization
peaks_reads <- read.table("~/Desktop/Rfiles/peak_unit/ZH11_RT_org.gff3", sep = "\t")
colnames(peaks_reads) <- c("chr","start","end","RT")
RT_freq <- peaks_reads %>% group_by(RT)%>%
  summarise(Sum = sum(end) - sum(start))%>%
  filter(!(RT %in% c("ESLS","ESMSLS","MSLS")))%>%
  mutate(percent = Sum/sum(Sum), RT = factor(RT,levels = c("ES", "ESMS", "MS","LS")))%>%
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
  filter(RT %in% c("ES","ESMS","MS","LS")) %>%
  arrange(chr, start, end) 

write.table(table(tf_location_nondup$TF,tf_location_nondup$RT),"~/Desktop/Rfiles/peak_unit/tf_location_classfied_org.bed",row.names = T,quote=F,sep = "\t",col.names = T)

tf_location_nondup <- as.data.frame(read.table("/Users/lzz/Desktop/Rfiles/peak_unit/tf_location_classfied_org.bed",header = T))
# real_bs <- read_xlsx("~/Desktop/Rfiles/peak_unit/x39_49_wx11_heat.xlsx")
# tf_location_nondup <- rbind(tf_location_nondup,real_bs)
tf_location_nondup <- tf_location_nondup[,c("ES","ESMS","MS","LS")]
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
data.1 <- as.data.frame(decostand(tf_location_nondup,"standardize",MARGIN = 1)) 

# Identify the replication phase where each TF has the highest binding value
max_binding_stage <- colnames(data.1)[apply(data.1, 1, which.max)]

# Convert to factor with predefined order
stage_order <- c("ES", "ESMS", "MS", "LS")
max_binding_stage <- factor(max_binding_stage, levels = stage_order)

# Sort transcription factors by their primary binding phase
sorted_indices <- order(max_binding_stage)

# Reorder data accordingly
data.1 <- data.1[sorted_indices, ]


pic_heatmap<-pheatmap(data.1,show_rownames = FALSE,show_colnames = TRUE,
                      display_numbers = matrix(ifelse(data.1 > 2, ""," "), 
                                               nrow(data.1)),cluster_rows = FALSE)
#,color = colorRampPalette(c("#547297", "#8C9EBA", "#D9E0E7","#F3DBD6","#DA8F87"
#                           ,"#D54846"))(100))
# ggsave("~/Documents/GitHub/repli-ATAC-seq/output/Figures/tf_location_ZH11_4RT.pdf", pic_heatmap , width = 8, height = 5)

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

# Identify the replication phase where each TF has the highest binding value
max_binding_stage <- colnames(a[,5:8])[apply(a[,5:8], 1, which.max)]

# Convert to factor with predefined order
stage_order <- c("E", "EM", "M", "L")
max_binding_stage <- factor(max_binding_stage, levels = stage_order)

# Sort transcription factors by their primary binding phase
sorted_indices <- order(max_binding_stage)

# Reorder data accordingly
a <- a[sorted_indices, ]
a$max_binding_stage <- colnames(a[,5:8])[apply(a[,5:8], 1, which.max)]
# write.table(a,"~/Documents/GitHub/repli-ATAC-seq/output/tf_location_NIP_4RT_all.csv",quote=F,sep = ",",row.names = F)

# Define the order of max_binding_stage
stage_order <- c("E", "EM", "M", "L")

# Define custom colors for max_binding_stage
stage_colors <- c("E" = "#2C5F9E", 
                  "EM" = "#68A0D8", 
                  "M" = "#95BE6C", 
                  "L" = "#E68364")

# Select the top 5 most frequent TF_family_TAIR for each max_binding_stage, excluding NA and "Unknown"
top_tf_families <- a %>%
  filter(!is.na(TF_family_TAIR), TF_family_TAIR != "Unknown") %>%  # Remove NA and "Unknown"
  group_by(max_binding_stage, TF_family_TAIR) %>%
  summarise(count = n(), .groups = "drop") %>%  # Count occurrences
  filter(count >= 5) %>%  # Exclude those with frequency less than 5
  arrange(max_binding_stage, desc(count)) %>%
  group_by(max_binding_stage) %>%
  slice_head(n = 5)  # Select the top 5 per stage

# Compute the mean values for E, EM, M, L
avg_binding_values <- a %>%
  group_by(max_binding_stage, TF_family_TAIR) %>%
  summarise(across(c(E, EM, M, L), mean, na.rm = TRUE), .groups = "drop")

# Keep only the top 5 TF_family_TAIR for each max_binding_stage
final_data <- top_tf_families %>%
  inner_join(avg_binding_values, by = c("max_binding_stage", "TF_family_TAIR"))

# Ensure max_binding_stage follows the specified order
final_data$max_binding_stage <- factor(final_data$max_binding_stage, levels = stage_order)

# Generate dotplot with custom colors and size reflecting count
p <- ggplot(final_data, aes(x = max_binding_stage, y = TF_family_TAIR)) +
  geom_point(aes(size = count, color = max_binding_stage), alpha = 0.8) +  # 透明度优化
  scale_color_manual(values = stage_colors) +  # 自定义颜色
  scale_size(range = c(3, 12)) +  # 调整大小范围
  theme_classic() +  # 使用更清晰的主题
  labs(
    title = "Top 5 TF families in Each stage",
    x = "Binding Stage",
    y = "TF Family",
    size = "Frequency",
    color = "Binding Stage"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_line(size = 0.8),  # 添加坐标轴
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # 添加柔和网格线
    legend.position = "right"
  ) +
  guides(size = guide_legend(order = 1), color = guide_legend(order = 2))  # 调整图例顺序
ggsave("~/Documents/GitHub/repli-ATAC-seq/output/Figures/Dotplot.pdf", p , width = 8, height = 5)
