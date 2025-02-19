library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/macs2/macs2_noctrl_p0.01/rm_noisepeak/results/")

# **TPM Normalization**
normalize_tpm_cr <- function(data) {
  data$RegionLength <- data$end - data$start + 1
  data[, 4:18] <- lapply(data[, 4:18], function(col) {
    rpk <- col / (data$RegionLength / 1000)
    scaling_factor <- sum(rpk) / 1e6
    tpm <- rpk / scaling_factor
    return(tpm)
  })
  return(data)
}

normalize_tpm_wt <- function(data) {
  data$RegionLength <- data$end - data$start + 1
  data[, 4:21] <- lapply(data[, 4:21], function(col) {
    rpk <- col / (data$RegionLength / 1000)
    scaling_factor <- sum(rpk) / 1e6
    tpm <- rpk / scaling_factor
    return(tpm)
  })
  return(data)
}

# **Load Data**
wt <- read.table("./peaks_re_quan_org.txt", header = FALSE)
colnames(wt) <- c("chr", "start", "end",
                  "ZH11-3-ES", "ZH11-3-MS", "ZH11-3-LS",
                  "ZH11-4-ES", "ZH11-4-MS", "ZH11-4-LS",
                  "ZH11-2-G1", "ZH11-2-ES", "ZH11-2-MS", "ZH11-2-LS",
                  "ZH11-1-G1", "ZH11-1-ES", "ZH11-1-MS", "ZH11-1-LS",
                  "NIP-1-G1", "NIP-1-ES", "NIP-1-MS", "NIP-1-LS")

CR <- read.table("./peaks_re_quan_CR_org.txt", header = FALSE)
colnames(CR) <- c("chr", "start", "end",
                  "xw11-OE-ES", "xw11-OE-MS", "xw11-OE-LS",
                  "xw11-CR-ES", "xw11-CR-MS", "xw11-CR-LS",
                  "x39-CR-ES", "x39-CR-MS", "x39-CR-LS",
                  "ZH11-3-ES", "ZH11-3-MS", "ZH11-3-LS",
                  "ZH11-4-ES", "ZH11-4-MS", "ZH11-4-LS")

# **TPM Normalization**
wt_tpm <- normalize_tpm_wt(wt)
CR_tpm <- normalize_tpm_cr(CR)

# **Load Peaks and CUT&Tag Data**
peaks_ES <- read.table("/storage/liuxiaodongLab/liaozizhuo/Projects/repli-ATAC-seq/reference/ZH11_RT_all_org.gff3", header = FALSE)
colnames(peaks_ES) <- c("chr", "start", "end", "feature")
peaks_ES <- peaks_ES %>% filter(feature == "ES") %>% select(chr, start, end)

cuttag <- read.table("/storage/liuxiaodongLab/liaozizhuo/Projects/cuttag/macs2/macs2_p1e-5/xw11_cut_out/xw11_RT_top1000_org.bed", header = FALSE)
colnames(cuttag) <- c("chr", "start", "end", paste0("extra", 1:10), "cuttag_label")
cuttag <- cuttag %>% filter(cuttag_label == "ES") %>% select(chr, start, end)

# **Merge Data**
common_peaks <- inner_join(wt_tpm, CR_tpm, by = c("chr", "start", "end")) %>%
  inner_join(peaks_ES, by = c("chr", "start", "end"))

# **Compute Mean WT TPM**
wt_avg <- common_peaks %>%
  mutate(WT_avg_ES_TPM = rowMeans(select(., contains("ES")), na.rm = TRUE)) %>%
  select(chr, start, end, WT_avg_ES_TPM)

# **Merge WT Average with CR/OE Data**
combined_data <- inner_join(wt_avg, common_peaks, by = c("chr", "start", "end")) %>%
  mutate(cuttag_overlap = pmap_lgl(list(chr, start, end), function(c, s, e) {
    any(cuttag$chr == c & cuttag$start <= e & cuttag$end >= s)
  }))

# **Convert to Long Format (Only ES)**
plot_data <- combined_data %>%
  pivot_longer(cols = c(WT_avg_ES_TPM, starts_with("xw11"), starts_with("x39")), 
               names_to = "Sample", values_to = "TPM") %>%
  filter(grepl("ES", Sample)) %>%
  mutate(Batch = case_when(
    grepl("x39-CR", Sample) ~ "x39-CR",
    grepl("xw11-CR", Sample) ~ "xw11-CR",
    grepl("xw11-OE", Sample) ~ "xw11-OE",
    grepl("WT_avg", Sample) ~ "WT",
    TRUE ~ "Other"
  ))

# **Peak ID Assignment**
plot_data <- plot_data %>%
  group_by(Batch) %>%
  mutate(Peak_ID = row_number()) %>%
  ungroup()

# **Color Settings**
num_samples <- length(unique(plot_data$Sample))
color_palette <- scales::hue_pal()(num_samples)

# **Main Plots**
p1 <- ggplot(plot_data, aes(x = Peak_ID, y = TPM, color = Sample, group = Sample)) +
  geom_line(linewidth = 0.7, alpha = 0.5) +
  facet_wrap(~ Batch, scales = "free_y") +
  scale_color_manual(values = color_palette) +
  labs(x = "Peak Regions", y = "TPM",
       title = "Early Stage (ES) - WT vs Mutant & Overexpression") +
  theme_minimal() + ylim(0, 300)

p2 <- ggplot(plot_data %>% filter(cuttag_overlap), aes(x = Peak_ID, y = TPM, color = Sample, group = Sample)) +
  geom_line(linewidth = 0.7, alpha = 0.5) +
  facet_wrap(~ Batch, scales = "free_y") +
  scale_color_manual(values = color_palette) +
  labs(x = "CutTag Peaks", y = "TPM",
       title = "Early Stage (ES) - CutTag Peaks") +
  theme_minimal() + ylim(0, 300)

# **Scatter Plots (Log Scale)**
tpm_min <- min(combined_data$WT_avg_ES_TPM, combined_data$`xw11-CR-ES`, combined_data$`xw11-OE-ES`, na.rm = TRUE) + 1e-3
tpm_max <- max(combined_data$WT_avg_ES_TPM, combined_data$`xw11-CR-ES`, combined_data$`xw11-OE-ES`, na.rm = TRUE) + 1e-3

# **Scatter Plot: WT vs CR (Log Scale)**
p3 <- ggplot(combined_data, aes(x = WT_avg_ES_TPM + 1e-3, y = `xw11-CR-ES` + 1e-3, color = cuttag_overlap)) +
  geom_point(alpha = ifelse(scatter_data$cuttag_overlap, 1, 0.2), size = ifelse(scatter_data$cuttag_overlap, 2, 1)) +
  scale_color_manual(values = c("gray60", "red")) +  
  scale_x_log10(limits = c(tpm_min, tpm_max)) +  
  scale_y_log10(limits = c(tpm_min, tpm_max)) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "WT ES TPM (log10)", y = "CR ES TPM (log10)", title = "WT vs. CR") +
  coord_fixed(ratio = 1) +  
  theme_minimal()

# **Scatter Plot: WT vs OE (Log Scale)**
p4 <- ggplot(combined_data, aes(x = WT_avg_ES_TPM + 1e-3, y = `xw11-OE-ES` + 1e-3, color = cuttag_overlap)) +
  geom_point(alpha = ifelse(scatter_data$cuttag_overlap, 1, 0.2), size = ifelse(scatter_data$cuttag_overlap, 2, 1)) +
  scale_color_manual(values = c("gray60", "red")) +  
  scale_x_log10(limits = c(tpm_min, tpm_max)) +  
  scale_y_log10(limits = c(tpm_min, tpm_max)) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "WT ES TPM (log10)", y = "OE ES TPM (log10)", title = "WT vs. OE") +
  coord_fixed(ratio = 1) +  
  theme_minimal()

print(p3)
print(p4)

