# Define four matrices
wt <- matrix(c(165, 75, 51, 181,
               12820, 18017, 29502, 28550), nrow = 2, byrow = TRUE)

xw11_OE <- matrix(c(172, 133, 55, 112,
                    13869, 20473, 31967, 22580), nrow = 2, byrow = TRUE)

xw11_CR <- matrix(c(153, 163, 58, 98,
                    15086, 24780, 33055, 15968), nrow = 2, byrow = TRUE)

x39_CR <- matrix(c(207, 103, 61, 101,
                   16442, 20882, 33870, 17695), nrow = 2, byrow = TRUE)

# Set column and row names
colnames(wt) <- colnames(xw11_OE) <- colnames(xw11_CR) <- colnames(x39_CR) <- c("ES", "MS", "LS", "Other")
rownames(wt) <- rownames(xw11_OE) <- rownames(xw11_CR) <- rownames(x39_CR) <- c("Target", "Background")

# Function to calculate enrichment and p-values (including "Other")
calculate_enrichment_pvalues <- function(matrix) {
  col_totals <- colSums(matrix)
  row_totals <- rowSums(matrix)
  grand_total <- sum(matrix)
  
  # Compute enrichment for all categories including "Other"
  enrichment <- matrix[1, ] / (col_totals * row_totals[1] / grand_total)
  
  # Compute p-values using Fisher's exact test for each category vs all others
  p_values <- sapply(1:4, function(i) {
    other_cols <- setdiff(1:4, i)  # Select all columns except the one being tested
    test_matrix <- matrix[, c(i, other_cols[1])]  # Compare category to the first "other" column
    fisher.test(test_matrix, simulate.p.value = TRUE)$p.value
  })
  
  # Store results in a data frame
  result <- data.frame(
    Category = colnames(matrix),
    Enrichment = enrichment,
    P_Value = p_values
  )
  
  return(result)
}

# Apply function to all matrices
wt_results <- calculate_enrichment_pvalues(wt)
xw11_OE_results <- calculate_enrichment_pvalues(xw11_OE)
xw11_CR_results <- calculate_enrichment_pvalues(xw11_CR)
x39_CR_results <- calculate_enrichment_pvalues(x39_CR)

# Print results
print("WT Results:")
print(wt_results)

print("xw11_OE Results:")
print(xw11_OE_results)

print("xw11_CR Results:")
print(xw11_CR_results)

print("x39_CR Results:")
print(x39_CR_results)