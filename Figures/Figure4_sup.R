# Load necessary libraries
library(dplyr)
library(pheatmap)

# Step 1: Read TF binding data
tf_family <- read.table("/Users/lzz/Desktop/Rfiles/fimo.bed",header = T ,sep = "\t")
tf_family <- tf_family[,c(2,3)]
tf_family <- distinct(tf_family,motif_alt_id ,TF_family,.keep_all= TRUE)
colnames(tf_family)[1] <- "TF"

tf_table <- read.table("~/Desktop/Rfiles/peak_unit/tf_location_classfied_org.bed",
                       header = TRUE, sep = "\t", row.names = 1)

# Keep only the four relevant columns: ES, ESMS, MS, LS
tf_table <- tf_table[, c("ES", "ESMS", "MS", "LS")]
colnames(tf_table) <- c("E", "EM", "M", "L")  # Rename columns for consistency

# Step 2: Read region length data for each stage (E, EM, M, L)
RT_length <- tibble::tibble(
  RT = c("E", "EM", "M", "L"),
  Sum = c(10008822, 1659300, 14546125, 18878938)
)
RT_lengths <- setNames(RT_length$Sum, RT_length$RT)

# Step 3: Compute expected values
# Expected values assume TF binding follows the proportion of the region length
tf_totals <- rowSums(tf_table)
expected_values <- as.data.frame(sapply(RT_lengths, function(x) x / sum(RT_lengths) * tf_totals))
rownames(expected_values) <- rownames(tf_table)

# Step 4: Compute Fold Change (with pseudo-count to prevent zero division)
epsilon <- 1
fold_change <- (tf_table + epsilon) / (expected_values + epsilon)

# Step 5: Compute log2 Fold Change as an enrichment score
enrichment_score <- log2(fold_change)

# Step 6: Compute Running Sum (GSEA method)
running_sum <- function(fc_values) {
  sorted_values <- sort(fc_values, decreasing = TRUE)  # Sort TFs by log2FC (descending)
  running_sum_values <- cumsum(sorted_values - mean(sorted_values))  # Cumulative sum
  return(running_sum_values)
}

# Apply Running Sum computation for each stage (E, EM, M, L)
running_sums <- apply(enrichment_score, 2, running_sum)

# Step 7: Perform Permutation Test to compute p-values
set.seed(123)  # Set seed for reproducibility
num_permutations <- 1000
permuted_scores <- matrix(NA, nrow = num_permutations, ncol = ncol(enrichment_score))

# Generate permuted datasets and compute random enrichment scores
for (i in 1:num_permutations) {
  permuted_tf <- sample(rownames(enrichment_score))  # Shuffle TF names randomly
  permuted_es <- enrichment_score[permuted_tf, ]  # Rearrange data accordingly
  permuted_scores[i, ] <- colMeans(permuted_es)  # Compute mean enrichment score
}

# Compute p-values by comparing observed vs. permuted enrichment scores
p_values <- apply(enrichment_score, 2, function(obs) {
  mean(abs(obs) >= abs(permuted_scores))
})

# Perform FDR correction (Benjamini-Hochberg method)
fdr_values <- p.adjust(p_values, method = "BH")

# Step 8: Identify the most enriched stage for each TF
max_enriched_stage <- apply(enrichment_score, 1, function(x) colnames(enrichment_score)[which.max(x)])

# Create a final table of enriched TFs
tf_enrichment_final <- data.frame(
  TF = rownames(enrichment_score),
  Max_Enriched_Stage = max_enriched_stage,
  Max_Enrichment = apply(enrichment_score, 1, max)
)

# Step 9: Sort TFs by enrichment stage (E → EM → M → L) and descending enrichment score
stage_order <- c("E", "EM", "M", "L")
tf_enrichment_final$Max_Enriched_Stage <- factor(tf_enrichment_final$Max_Enriched_Stage, levels = stage_order)
tf_enrichment_final <- tf_enrichment_final %>% arrange(Max_Enriched_Stage, desc(Max_Enrichment))

# Step 10: Reorder the enrichment score matrix based on TF ranking
enrichment_score_sorted <- enrichment_score[match(tf_enrichment_final$TF, rownames(enrichment_score)), ]

# Step 11: Generate heatmap for TF-level enrichment scores (without TF names)
pheatmap(enrichment_score_sorted, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "TF-Level Enrichment Score",  # Heatmap title
         cellwidth = 90,  # Set cell width (adjust as needed)
         cellheight = 0.6,  # Set cell height (adjust as needed)
         show_rownames = FALSE,  # Hide TF names
         legend_labels = "Enrichment Score")  # Manually setting legend title

# Step 12: Compute enrichment scores for TF families
# Load TF-to-family mapping file (TF family information)
# Example file format: first column (TF), second column (TF family)

# Step 13: Compute enrichment scores at the TF family level
# Merge TF binding data with TF family mapping
# Compute TF family-level enrichment scores
tf_family_enrichment <- enrichment_score %>%
  as.data.frame() %>%
  rownames_to_column("TF") %>%
  inner_join(tf_family, by = "TF") %>%  # Map TFs to their families
  group_by(TF_family) %>%
  summarise(across(E:L, mean, na.rm = TRUE)) %>%  # Compute average enrichment score per family
  column_to_rownames("TF_family")  # Set TF_family as row names

# Ensure column order is E → EM → M → L
tf_family_enrichment <- tf_family_enrichment[, c("E", "EM", "M", "L")]

# Identify the most enriched stage for each TF family
max_enriched_stage <- apply(tf_family_enrichment, 1, function(x) colnames(tf_family_enrichment)[which.max(x)])

# Create a sorted table
tf_family_sorted <- data.frame(
  TF_Family = rownames(tf_family_enrichment),
  Max_Enriched_Stage = max_enriched_stage,
  Max_Enrichment = apply(tf_family_enrichment, 1, max)
)

# Define the order of stages
stage_order <- c("E", "EM", "M", "L")

# Convert the Max_Enriched_Stage column to a factor with the defined order
tf_family_sorted$Max_Enriched_Stage <- factor(tf_family_sorted$Max_Enriched_Stage, levels = stage_order)

# Sort TF families first by stage, then by highest enrichment score
tf_family_sorted <- tf_family_sorted %>% arrange(Max_Enriched_Stage, desc(Max_Enrichment))

# Reorder the enrichment score matrix based on the sorted TF families
tf_family_enrichment_sorted <- tf_family_enrichment[match(tf_family_sorted$TF_Family, rownames(tf_family_enrichment)), ]

# Generate heatmap with sorted TF family enrichment scores
pheatmap(as.matrix(tf_family_enrichment_sorted), 
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "TF Family-Level Enrichment Score",  # Heatmap title
         cellwidth = 40,  # Set cell width (adjust as needed)
         cellheight = 15,  # Set cell height (adjust as needed)
         show_rownames = TRUE)  # Show family names
