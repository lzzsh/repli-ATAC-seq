#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.umask("0077")
suppressPackageStartupMessages(library(edgeR))

task_root <- "/storage2/liuxiaodongLab/liaozizhuo/Projects/extended/2026-08-26_009_cpp8-allele-path-replication"
counts_file <- "/storage2/liuxiaodongLab/liaozizhuo/Projects/RNA-seq/star_all_rawdata/featureCounts_out/gene_counts.txt"
out_file <- file.path(task_root, "output", "rna_joint_edger_all_genes_009.tsv.gz")
sample_file <- file.path(task_root, "data", "rna_sample_manifest_009.tsv")
norm_file <- file.path(task_root, "output", "rna_tmm_factors_009.tsv")
log_file <- file.path(task_root, "output", "logs", "02_build_joint_rna_contrasts.log")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(sample_file), recursive = TRUE, showWarnings = FALSE)
if (!startsWith(normalizePath(dirname(out_file)), normalizePath(task_root))) {
  stop("Refusing output outside task root")
}

raw <- read.delim(counts_file, header = TRUE, comment.char = "#", check.names = FALSE)
count_matrix <- as.matrix(raw[, 7:ncol(raw)])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- raw$Geneid

original_names <- c(
  "TCX2-3-KO.1", "WT.1", "WT.2", "SOL1-5-KO.1", "SOL1-5-KO.2",
  "SOL1-8-KO.1", "SOL1-8-KO.2", "TCX2-1-KO.1", "TCX2-1-KO.2", "TCX2-3-KO.2"
)
if (ncol(count_matrix) != length(original_names)) stop("Unexpected RNA count column number")
colnames(count_matrix) <- original_names

selected <- c("WT.1", "WT.2", "TCX2-1-KO.1", "TCX2-1-KO.2", "TCX2-3-KO.1", "TCX2-3-KO.2")
groups <- factor(c("WT", "WT", "CPP8_1", "CPP8_1", "CPP8_3", "CPP8_3"), levels = c("WT", "CPP8_1", "CPP8_3"))
count_sub <- count_matrix[, selected, drop = FALSE]

y0 <- DGEList(counts = count_sub, group = groups)
# A single three-group filter retains genes expressed in only one mutant and is
# too permissive for a shared two-contrast path universe.  Require eligibility
# in both pairwise WT-versus-allele comparisons, then estimate both contrasts
# together in the same six-sample model.
idx1 <- groups %in% c("WT", "CPP8_1")
idx3 <- groups %in% c("WT", "CPP8_3")
y_filter1 <- DGEList(counts = count_sub[, idx1, drop = FALSE], group = droplevels(groups[idx1]))
y_filter3 <- DGEList(counts = count_sub[, idx3, drop = FALSE], group = droplevels(groups[idx3]))
keep1 <- filterByExpr(y_filter1, group = droplevels(groups[idx1]))
keep3 <- filterByExpr(y_filter3, group = droplevels(groups[idx3]))
keep <- keep1 & keep3
y <- y0[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

contrast_1 <- makeContrasts(CPP8_1 - WT, levels = design)
contrast_3 <- makeContrasts(CPP8_3 - WT, levels = design)
tab1 <- topTags(glmQLFTest(fit, contrast = contrast_1), n = Inf, sort.by = "none")$table
tab3 <- topTags(glmQLFTest(fit, contrast = contrast_3), n = Inf, sort.by = "none")$table

norm_cpm <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
wt_mean_cpm <- rowMeans(norm_cpm[, c("WT.1", "WT.2"), drop = FALSE])

common <- data.frame(
  gene_id = rownames(y$counts),
  wt_mean_cpm = wt_mean_cpm[rownames(y$counts)],
  cpp8_1_logFC = tab1[rownames(y$counts), "logFC"],
  cpp8_1_logCPM = tab1[rownames(y$counts), "logCPM"],
  cpp8_1_PValue = tab1[rownames(y$counts), "PValue"],
  cpp8_1_FDR = tab1[rownames(y$counts), "FDR"],
  cpp8_3_logFC = tab3[rownames(y$counts), "logFC"],
  cpp8_3_logCPM = tab3[rownames(y$counts), "logCPM"],
  cpp8_3_PValue = tab3[rownames(y$counts), "PValue"],
  cpp8_3_FDR = tab3[rownames(y$counts), "FDR"],
  row.names = NULL,
  check.names = FALSE
)

con <- gzfile(out_file, "wt")
write.table(common, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

sample_manifest <- data.frame(
  sample = selected,
  condition = as.character(groups),
  source_count_column = match(selected, original_names),
  featureCounts_bam_column = colnames(raw)[6 + match(selected, original_names)]
)
write.table(sample_manifest, sample_file, sep = "\t", quote = FALSE, row.names = FALSE)

norm <- data.frame(
  sample = colnames(y),
  condition = as.character(groups),
  library_size = y$samples$lib.size,
  tmm_norm_factor = y$samples$norm.factors,
  effective_library_size = y$samples$lib.size * y$samples$norm.factors
)
write.table(norm, norm_file, sep = "\t", quote = FALSE, row.names = FALSE)

log_lines <- c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("R=", R.version.string),
  paste0("edgeR=", as.character(packageVersion("edgeR"))),
  paste0("counts_file=", counts_file),
  paste0("genes_raw=", nrow(count_matrix)),
  paste0("genes_pass_cpp8_1_pairwise_filter=", sum(keep1)),
  paste0("genes_pass_cpp8_3_pairwise_filter=", sum(keep3)),
  paste0("genes_joint_filter=", nrow(common)),
  paste0("design_columns=", paste(colnames(design), collapse = ",")),
  paste0("output=", out_file)
)
writeLines(log_lines, log_file)
Sys.chmod(c(out_file, sample_file, norm_file, log_file), mode = "0600")

cat("Joint filtered genes:", nrow(common), "\n")
cat("Wrote", out_file, "\n")
