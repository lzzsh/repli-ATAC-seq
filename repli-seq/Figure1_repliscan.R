# ========== 加载依赖包 ==========
library(dplyr)
library(ggplot2)
library(ggrepel)

# ========== 定义绘图函数 ==========
plot_rt_from_gff <- function(gff_path, sample_name, save_dir = "./") {
  gff <- read.table(gff_path, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
  gff <- gff[, c(1, 4, 5, 9)]
  colnames(gff) <- c("chr", "start", "end", "RT")
  
  gff$RT <- sapply(gff$RT, function(x) {
    tag <- unlist(strsplit(x, ";"))[2]
    sub(".*=", "", tag)
  })
  
  gff <- gff[!(gff$chr %in% c("chrUn", "chrSy")), ]
  gff <- gff[!is.na(gff$RT), ]
  gff$RT <- gsub("S", "", gff$RT)
  gff$length <- gff$end - gff$start + 1
  
  # --- 饼图 Figure 1A ---
  gff$RT <- factor(gff$RT, levels = c("EML", "EL", "L", "ML", "M", "EM", "E"))
  Figure1A <- ggplot(gff, aes(x = "", fill = RT)) +
    geom_bar(width = 1, aes(weight = length), stat = "count") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c(
      E = "#2C5F9E", EM = "#68A0D8", M = "#95BE6C",
      ML = "#E4B660", L = "#E68364", EL = "#B784A7", EML = "#9B7EB3"
    )) +
    labs(title = paste0(sample_name, " - Partition of genome replication")) +
    theme_void()
  
  ggsave(
    filename = file.path(save_dir, paste0(sample_name, "_RT_Pie.pdf")),
    plot = Figure1A, width = 8, height = 5
  )
  
  # --- 箱线图 Figure 1B ---
  gff$log2length <- log2(gff$length)
  gff$RT <- factor(gff$RT, levels = c("E", "EM", "M", "ML", "L", "EL", "EML"))
  Figure1B <- ggplot(data = gff, aes(x = RT, y = log2length)) +
    geom_boxplot(aes(color = RT), size = 1.1) +
    scale_color_manual(values = c(
      E = "#2C5F9E", EM = "#68A0D8", M = "#95BE6C",
      ML = "#E4B660", L = "#E68364", EL = "#B784A7", EML = "#9B7EB3"
    )) +
    theme_classic() +
    ggtitle(paste0(sample_name, " - RT segment size distribution")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("log2(Length)")
  
  ggsave(
    filename = file.path(save_dir, paste0(sample_name, "_RT_Boxplot.pdf")),
    plot = Figure1B, width = 8, height = 5
  )
}

# ========== 指定样本路径与名称 ==========
samples <- list(
  wt_norep = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_wt_norep.gff3",
  wt       = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_wt.gff3",
  sol1_5   = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_cr_sol1_5.gff3",
  sol1_8   = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_cr_sol1_8.gff3",
  tcx2_1   = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_cr_tcx2_1.gff3",
  tcx2_3   = "~/Desktop/repli-ATAC-seq/segamentation/ratio_segmentation_cr_tcx2_3.gff3"
)

# ========== 设置输出路径 ==========
output_dir <- "~/Documents/GitHub/repli-ATAC-seq/output/Figures"

# ========== 批量绘图 ==========
for (sample_name in names(samples)) {
  cat("Processing:", sample_name, "\n")
  plot_rt_from_gff(gff_path = samples[[sample_name]],
                   sample_name = sample_name,
                   save_dir = output_dir)
}