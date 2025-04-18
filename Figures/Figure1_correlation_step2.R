library(tidyverse)

# Step 1: List all raw count files from bedtools coverage
raw_files <- list.files("coverage_raw", pattern = "_raw.txt$", full.names = TRUE)

# Step 2: Read each file and assign column names
raw_list <- lapply(raw_files, function(file) {
  df <- read.table(file, header = FALSE, col.names = c("window", "count"))
  sample_name <- sub("_raw.txt", "", basename(file))  # Extract sample name
  colnames(df)[2] <- sample_name
  return(df)
})

# Step 3: Merge all sample data by genomic window
raw_merged <- Reduce(function(x, y) full_join(x, y, by = "window"), raw_list)

# Step 4: Replace missing values (NA) with 0 (no reads)
raw_merged[is.na(raw_merged)] <- 0

# Step 5: Extract count matrix and set row names as genomic windows
count_matrix <- raw_merged[,-1]  # Drop the 'window' column
rownames(count_matrix) <- raw_merged$window

# Step 6: Compute total read count (library size) for each sample
total_counts <- colSums(count_matrix)

# Step 7: Normalize to RPM (Reads Per Million)
rpm_matrix <- sweep(count_matrix, 2, total_counts / 1e6, FUN = "/")

# Step 8: Log-transform the RPM matrix: log2(RPM + 1)
log2rpm_matrix <- log2(rpm_matrix + 1)

# Step 9: Combine windows and transformed data, write to output
output_matrix <- data.frame(window = rownames(log2rpm_matrix), log2rpm_matrix)
write.table(output_matrix, file = "ATAC_log2RPM_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)