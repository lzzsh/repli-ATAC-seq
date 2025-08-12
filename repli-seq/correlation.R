library(tidyverse)
library(pheatmap)
library(ggrepel)

# ========== Step 1: Load data ==========
raw_matrix <- read.table("/Users/lzz/Desktop/repli-ATAC-seq/repli-seq-CR/raw_counts_repli_cr.txt", sep = "\t")

# ========== Step 2: Prepare count matrix ==========
count_matrix <- raw_matrix[,-c(1,2,3)]
rownames(count_matrix) <- raw_matrix$window
colnames(count_matrix) <- c(
  "WT-1-G1", "WT-2-G1", "WT-1-ES", "WT-2-ES", "WT-1-MS", "WT-1-LS",
  "sol1_5-1-G1", "sol1_5-2-G1", "sol1_5-1-ES", "sol1_5-2-ES", "sol1_5-1-MS", "sol1_5-1-LS",
  "sol1_8-1-G1", "sol1_8-2-G1", "sol1_8-1-ES", "sol1_8-2-ES", "sol1_8-1-MS", "sol1_8-1-LS",
  "tcx2_1-1-G1", "tcx2_1-2-G1", "tcx2_1-1-ES", "tcx2_1-2-ES", "tcx2_1-1-MS", "tcx2_1-1-LS",
  "tcx2_3-1-G1", "tcx2_3-1-ES", "tcx2_3-1-MS", "tcx2_3-1-LS"
)

# ========== Step 3: Normalize to RPM ==========
total_counts <- colSums(count_matrix)
rpm_matrix <- sweep(count_matrix, 2, total_counts / 1e6, FUN = "/")

# ========== Step 4: Log2(RPM + 1) ==========
log2rpm_matrix <- log2(rpm_matrix + 1) + 1e-6  # Adding small constant for numerical stability

# ========== Step 5: Normalize each ES/MS/LS by corresponding G1 ==========
log2rpm_normalized <- log2rpm_matrix  # Copy matrix
# 
# for (i in 1:ncol(log2rpm_matrix)) {
#   sample_name <- colnames(log2rpm_matrix)[i]
#   
#   # Skip G1 samples
#   if (grepl("-G1$", sample_name)) next
#   
#   # Try to find the matching G1 sample by replacing suffix
#   g1_name <- sub("-(ES|MS|LS)$", "-G1", sample_name)
#   
#   if (g1_name %in% colnames(log2rpm_matrix)) {
#     log2rpm_normalized[, i] <- log2rpm_matrix[, i] - log2rpm_matrix[, g1_name]
#   } else {
#     warning(paste("No matching G1 for:", sample_name))
#     log2rpm_normalized[, i] <- NA
#   }
# }

# ========== Step 6: Compute sample correlation ==========
cor_matrix <- cor(log2rpm_normalized, use = "pairwise.complete.obs", method = "pearson")

# ========== Step 7: Correlation heatmap ==========
my_color <- colorRampPalette(c("red", "#FFF000", "white"))(100)
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

# ========== Step 8: PCA ==========
# Remove constant rows
row_sd <- apply(log2rpm_normalized, 1, sd, na.rm = TRUE)
log2rpm_filtered <- log2rpm_normalized[row_sd > 0, ]

pca_res <- prcomp(t(log2rpm_filtered), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
pca_df$stage <- sub(".*-(ES|MS|LS|G1)$", "\\1", pca_df$Sample)
time_colors <- c("ES" = "red", "MS" = "gold", "LS" = "blue", "G1" = "grey")
explained_var <- summary(pca_res)$importance[2, 1:2] * 100

pdf("/Users/lzz/Documents/GitHub/repli-ATAC-seq/output/Figures/ATAC_log2FC_pca.pdf", width = 8, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, color = stage, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values = time_colors) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "%)")
  ) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
