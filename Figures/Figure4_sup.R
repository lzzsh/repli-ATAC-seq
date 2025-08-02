# Load necessary libraries
library(dplyr)
library(pheatmap)

# Step 1: Read TF binding data
tf_family <- read.table("/Users/lzz/Desktop/Rfiles/fimo.bed", header = T, sep = "\t")
tf_family <- tf_family[, c(2, 3)]
tf_family <- distinct(tf_family, motif_alt_id, TF_family, .keep_all = TRUE)
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
tf_totals <- rowSums(tf_table)
expected_values <- as.data.frame(sapply(RT_lengths, function(x) x / sum(RT_lengths) * tf_totals))
rownames(expected_values) <- rownames(tf_table)

# Step 4: Compute Fold Change and enrichment score
epsilon <- 1
fold_change <- (tf_table + epsilon) / (expected_values + epsilon)
enrichment_score <- log2(fold_change)

# Step 5: Analyze individual TFs
max_enriched_stage <- apply(enrichment_score, 1, function(x) colnames(enrichment_score)[which.max(x)])

# Step 6: Sort TFs by enrichment stage and descending enrichment score
stage_order <- c("E", "EM", "M", "L")

# Create comprehensive TF results (combining results and matrix)
tf_results_comprehensive <- data.frame(
  TF = rownames(enrichment_score),
  Max_Enriched_Stage = factor(max_enriched_stage, levels = stage_order),
  Max_Enrichment = apply(enrichment_score, 1, max),
  E_Score = enrichment_score[, "E"],
  EM_Score = enrichment_score[, "EM"], 
  M_Score = enrichment_score[, "M"],
  L_Score = enrichment_score[, "L"],
  stringsAsFactors = FALSE
) %>% arrange(Max_Enriched_Stage, desc(Max_Enrichment))

# Step 7: Generate heatmap for TF-level enrichment scores
enrichment_score_sorted <- enrichment_score[match(tf_results_comprehensive$TF, rownames(enrichment_score)), ]

pheatmap(enrichment_score_sorted, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "TF-Level Enrichment Score",
         cellwidth = 90, cellheight = 0.6,
         show_rownames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100))

# Step 8: Compute enrichment scores at the TF family level
tf_family_enrichment <- enrichment_score %>%
  as.data.frame() %>%
  rownames_to_column("TF") %>%
  inner_join(tf_family, by = "TF") %>%
  group_by(TF_family) %>%
  summarise(across(E:L, mean, na.rm = TRUE), 
            TF_count = n(), .groups = 'drop') %>%
  column_to_rownames("TF_family")

# Ensure column order is E → EM → M → L
tf_family_enrichment_scores <- tf_family_enrichment[, c("E", "EM", "M", "L")]
tf_family_counts <- tf_family_enrichment[, "TF_count"]

# Identify the most enriched stage for each TF family
max_enriched_stage_family <- apply(tf_family_enrichment_scores, 1, function(x) colnames(tf_family_enrichment_scores)[which.max(x)])

# Create comprehensive TF family results (combining results and matrix)
tf_family_results_comprehensive <- data.frame(
  TF_Family = rownames(tf_family_enrichment_scores),
  Max_Enriched_Stage = factor(max_enriched_stage_family, levels = stage_order),
  Max_Enrichment = apply(tf_family_enrichment_scores, 1, max),
  TF_Count = tf_family_counts,
  E_Score = tf_family_enrichment_scores[, "E"],
  EM_Score = tf_family_enrichment_scores[, "EM"],
  M_Score = tf_family_enrichment_scores[, "M"],
  L_Score = tf_family_enrichment_scores[, "L"],
  stringsAsFactors = FALSE
) %>% arrange(Max_Enriched_Stage, desc(Max_Enrichment))

# Step 9: Generate heatmap with sorted TF family enrichment scores
tf_family_enrichment_sorted <- tf_family_enrichment_scores[match(tf_family_results_comprehensive$TF_Family, rownames(tf_family_enrichment_scores)), ]

pheatmap(t(as.matrix(tf_family_enrichment_sorted)), 
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "TF Family-Level Enrichment Score",
         cellwidth = 15, cellheight = 40,
         angle_col = 90,
         show_rownames = TRUE)

# Step 10: Print summary results
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total TFs analyzed:", nrow(enrichment_score), "\n")
cat("Total TF families analyzed:", nrow(tf_family_enrichment_scores), "\n\n")


cat("\nTop 5 most enriched TFs per stage:\n")
for (stage in stage_order) {
  top5 <- tf_results_comprehensive %>% 
    filter(Max_Enriched_Stage == stage) %>% 
    head(5)
  cat(sprintf("\n%s stage:\n", stage))
  for (i in 1:min(5, nrow(top5))) {
    cat(sprintf("  %d. %s (%.3f)\n", i, top5$TF[i], top5$Max_Enrichment[i]))
  }
}

# Step 11: Save comprehensive results (2 files only)
write.csv(tf_results_comprehensive, "TF_enrichment_comprehensive.csv", row.names = FALSE)
write.csv(tf_family_results_comprehensive, "TF_family_enrichment_comprehensive.csv", row.names = FALSE)
