# Load required libraries
library(ChIPseeker)
library(ggplot2)
library(ggimage)
library(ggupset)
library(GenomicFeatures)
library(ggpubr)

# Set ChIPseeker options
options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)

# Load narrowPeak files
files <- list.files("~/Desktop/Rfiles/peaks/")
setwd("~/Desktop/Rfiles/peaks")

# Build TxDb from rice GFF3 annotation
ricegff <- makeTxDbFromGFF("~/Desktop/Rfiles/all_DIY.gff3")

# Annotate three selected samples (ZH11 group: ES, LS, MS)
peakAnno1 <- annotatePeak(files[[7]], TxDb = ricegff, tssRegion = c(-3000, 3000),
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))

peakAnno2 <- annotatePeak(files[[9]], TxDb = ricegff, tssRegion = c(-3000, 3000),
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))

peakAnno3 <- annotatePeak(files[[8]], TxDb = ricegff, tssRegion = c(-3000, 3000),
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))

# Standardize annotation names and reorder levels
rename_annotation <- function(df) {
  # Keep full promoter annotation; strip gene details from others
  df$annotation <- ifelse(grepl("^Promoter", df$annotation),
                          df$annotation,
                          gsub(" \\(.*", "", df$annotation))
  
  # Remap remaining categories to clean names
  annotation_map <- c(
    "5' UTR"             = "5' UTR",
    "3' UTR"             = "3' UTR",
    "Exon"               = "Exon",
    "Intron"             = "Intron",
    "Downstream"         = "Downstream",
    "Downstream (<=300bp)" = "Downstream",
    "Distal Intergenic"  = "Intergenic"
  )
  df$annotation <- recode(df$annotation, !!!annotation_map, .default = df$annotation)
  
  # Specify the desired order of categories
  levels_order <- c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)",
                    "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Intergenic")
  
  df$annotation <- factor(df$annotation, levels = levels_order)
  return(df)
}

# Define custom cold-tone color palette (with yellow highlight)
my_colors <- c(
  "Promoter (<=1kb)" = "#0B21AB",  # highlighted in deep blue
  "Promoter (1-2kb)" = "#5555FC",
  "Promoter (2-3kb)" = "#FFF8DD",  # soft yellow
  "5' UTR"           = "#FED127",
  "3' UTR"           = "#841885",
  "Exon"             = "#B2DF8A",
  "Intron"           = "#E0E0E0",
  "Downstream"       = "#999999",
  "Intergenic"       = "#CAB2D6"
)

# General function to create a pie chart
make_pie_plot <- function(peakAnnoObj, title_label) {
  df <- data.frame(annotation = peakAnnoObj@anno$annotation) %>%
    count(annotation) %>%
    mutate(percentage = n / sum(n) * 100)
  df <- rename_annotation(df)
  
  ggplot(df, aes(x = "", y = percentage, fill = annotation)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = my_colors) +
    theme_void() +
    ggtitle(title_label) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "none")
}

# Create pie charts for the 3 stages
p1 <- make_pie_plot(peakAnno1, "ES")
p2 <- make_pie_plot(peakAnno2, "MS")
p3 <- make_pie_plot(peakAnno3, "LS")

# Arrange all 3 pie charts into one figure with shared legend
fig <- ggarrange(p1, p2, p3,
                 ncol = 3, nrow = 1,
                 common.legend = TRUE, legend = "right")

# Display figure
print(fig)

# Save figure to PDF
ggsave("/Users/lzz/Documents/GitHub/repli-ATAC-seq/output/Figures/proportion_of_peaks.pdf", fig, width = 12, height = 5)


# Step 1: Compute relative position in gene body
calculate_relative_position <- function(peakAnno, label) {
  df <- as.data.frame(peakAnno@anno)
  df <- df %>% filter(distanceToTSS > 3000)
  df$rel_pos <- (df$start - df$geneStart) / (df$geneEnd - df$geneStart)
  df$rel_pos <- ifelse(df$strand == "-", 1 - df$rel_pos, df$rel_pos)
  df <- df %>% filter(rel_pos >= 0 & rel_pos <= 1)
  df$type <- label
  return(df)
}

# Apply to each sample
rel_ES <- calculate_relative_position(peakAnno1, "ES")
rel_MS <- calculate_relative_position(peakAnno2, "MS")
rel_LS <- calculate_relative_position(peakAnno3, "LS")

# Combine all
all_peaks_zh11 <- bind_rows(rel_ES, rel_MS, rel_LS)

# Step 2: Assign to gene body quintile bins
all_peaks_zh11 <- all_peaks_zh11 %>%
  mutate(bin = cut(rel_pos, breaks = seq(0, 1, by = 0.2),
                   labels = paste0("Q", 1:5), include.lowest = TRUE))

# Step 3: Count and normalize to percentage
plot_data_zh11 <- all_peaks_zh11 %>%
  count(type, bin) %>%
  group_by(type) %>%
  mutate(percentage = n / sum(n) * 100)

# Step 4: Plot line chart (5 points per sample)
type_colors_zh11 <- c("ES" = "#0B21AB", "MS" = "#33a02c", "LS" = "#6a3d9a")
plot_data_zh11$type <- factor(plot_data_zh11$type, levels = c("ES", "MS", "LS"))
p5d <- ggplot(plot_data_zh11, aes(x = bin, y = percentage, group = type, color = type)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = type_colors_zh11) +
  labs(x = "Gene Body Quintile", y = "Percent of Peaks",
       title = "Peak Distribution Across Gene Body (Q1–Q5)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())


# Show figure
print(p5d)

# Save to PDF
ggsave("/Users/lzz/Documents/GitHub/repli-ATAC-seq/output/Figures/Figure5D_ZH11_gene_body_distribution.pdf",
       p5d, width = 7, height = 6)
