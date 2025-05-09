library(tidyverse)
library(pheatmap)
library(ggrepel)

# Step 1: Set working directory and load raw count table
raw_matrix <- read.table("/Users/lzz/Desktop/repli-ATAC-seq/reference_genome/raw_counts.txt", sep = "\t") 

# Step 2: Extract sample matrix (drop "window" column)
count_matrix <- raw_matrix[,-c(1,2,3)]
rownames(count_matrix) <- raw_matrix$window
colnames(count_matrix) <- c(
  "ZH11-4-ES", "ZH11-4-MS", "ZH11-4-LS",
  "ZH11-3-ES", "ZH11-3-MS", "ZH11-3-LS",
  "ZH11-2-ES", "ZH11-2-MS", "ZH11-2-LS",
  "ZH11-1-ES", "ZH11-1-MS", "ZH11-1-LS",
  "NIP-1-MS", "NIP-1-LS", "NIP-1-ES")

# Step 3: Normalize to RPM (Reads Per Million)
total_counts <- colSums(count_matrix)
rpm_matrix <- sweep(count_matrix, 2, total_counts / 1e6, FUN = "/")

# Step 4: Log-transform the RPM matrix: log2(RPM + 1)
log2rpm_matrix <- log2(rpm_matrix + 1)

# Step 5: Compute Pearson correlation
cor_matrix <- cor(log2rpm_matrix, method = "pearson")

# Step 6: Save correlation matrix
write.table(cor_matrix, file = "ATAC_log2RPM_sample_correlation.txt", sep = "\t", quote = FALSE)

# Step 7: Visualize as heatmap
my_color <- colorRampPalette(c("red", "#FFF000", "white"))(100)
pdf("/Users/lzz/Documents/GitHub/repli-ATAC-seq/output/Figures/ATAC_sample_correlation_heatmap.pdf", width = 8, height = 6)
pheatmap(cor_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = my_color,
         cellwidth = 20,        
         cellheight = 20, 
         border_color = NA,           
         display_numbers = FALSE,
         fontsize_number = 10)
dev.off()

# Step 8: Filter out constant rows (zero variance features)
log2rpm_matrix_filtered <- log2rpm_matrix[apply(log2rpm_matrix, 1, function(x) sd(x) > 0), ]

# Step 9: Perform PCA on transposed filtered matrix
pca_res <- prcomp(t(log2rpm_matrix_filtered), scale. = TRUE)

# Step 10: Format PCA result
pca_df <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
pca_df$stage <- sub(".*-(ES|MS|LS)$", "\\1", pca_df$Sample)

# Step 11: Draw PCA plot
time_colors <- c("ES" = "red", "MS" = "gold", "LS" = "blue")

explained_var <- summary(pca_res)$importance[2, 1:2] * 100 

# Step 13: Build PCA plot with labeled axes
pdf("/Users/lzz/Documents/GitHub/repli-ATAC-seq/output/Figures/ATAC_sample_correlation_pca.pdf", width = 8, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, color = stage, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values = time_colors) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "%)"),
  ) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
