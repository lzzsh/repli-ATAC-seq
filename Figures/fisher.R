# fisher test

# x39
fimo_site <- read.table("/Users/lzz/Desktop/Rfiles/peak_unit/tf_motif.bed")

fimo_site_x39 <- fimo_site %>%
  filter(V7 == "HB-5") %>%
  distinct(V2,V3,.keep_all = TRUE) %>%
  setNames(c("chr","start","end")) %>%
  dplyr::select(1,2,3) %>%
  mutate(fimo = 1)

cut_site_x39 <- read.table("/Users/lzz/Desktop/Rfiles/cuttag/x39_RT.bed")
cut_site_x39 <- cut_site_x39 %>%
  distinct(V5,V6,.keep_all = TRUE) %>%
  dplyr::select(1,5,6) %>%
  setNames(c("chr","start","end")) %>%
  mutate(cut = 1)

peaks_site <- peaks_reads[,c(1,2,3)]
peaks_site <- merge(peaks_site,fimo_site_x39,by=c("chr","start","end"), all.x = TRUE)
peaks_site <- merge(peaks_site,cut_site_x39,by=c("chr","start","end"), all.x = TRUE)

n11 <- length(which(peaks_site$fimo & peaks_site$cut))
n10 <- length(which(peaks_site$fimo & is.na(peaks_site$cut)))
n01 <- length(which(is.na(peaks_site$fimo) & peaks_site$cut))
n00 <- length(which(is.na(peaks_site$fimo) & is.na(peaks_site$cut)))
fisher_result <- fisher.test(matrix(c(n11, n10, n01, n00), ncol = 2))
fisher_result
matrix(c(n11, n10, n01, n00), ncol = 2)

# x49
fimo_site <- read.table("/Users/lzz/Desktop/Rfiles/peak_unit/tf_motif.bed")

fimo_site_x49 <- fimo_site %>%
  filter(V7 == "HAT5") %>%
  distinct(start,end,.keep_all = TRUE) %>%
  setNames(c("chr","start","end")) %>%
  dplyr::select(1,2,3) %>%
  mutate(fimo = 1)

cut_site_x49 <- read.table("/Users/lzz/Desktop/Rfiles/cuttag/x49_RT.bed")
cut_site_x49 <- cut_site_x49 %>%
  distinct(V5,V6,.keep_all = TRUE) %>%
  dplyr::select(1,5,6) %>%
  setNames(c("chr","start","end")) %>%
  mutate(cut = 1)

peaks_site <- peaks_reads[,c(1,2,3)]
peaks_site <- merge(peaks_site,fimo_site_x49,by=c("chr","start","end"), all.x = TRUE)
peaks_site <- merge(peaks_site,cut_site_x49,by=c("chr","start","end"), all.x = TRUE)

n11 <- length(which(peaks_site$fimo & peaks_site$cut))
n10 <- length(which(peaks_site$fimo & is.na(peaks_site$cut)))
n01 <- length(which(is.na(peaks_site$fimo) & peaks_site$cut))
n00 <- length(which(is.na(peaks_site$fimo) & is.na(peaks_site$cut)))
fisher_result <- fisher.test(matrix(c(n11, n10, n01, n00), ncol = 2))
fisher_result
matrix(c(n11, n10, n01, n00), ncol = 2)

# xw11
fimo_site <- read.table("/Users/lzz/Desktop/Rfiles/peak_unit/tf_motif.bed")

fimo_site_xw11 <- fimo_site %>%
  filter(V7 == "WOX11") %>%
  distinct(V2,V3,.keep_all = TRUE) %>%
  setNames(c("chr","start","end")) %>%
  dplyr::select(1,2,3) %>%
  mutate(fimo = 1)

cut_site_xw11 <- read.table("/Users/lzz/Desktop/Rfiles/cuttag/xw11_RT.bed")
cut_site_xw11 <- cut_site_xw11 %>%
  distinct(V5,V6,.keep_all = TRUE) %>%
  dplyr::select(1,5,6) %>%
  setNames(c("chr","start","end")) %>%
  mutate(cut = 1)

peaks_site <- peaks_reads[,c(1,2,3)]
peaks_site <- merge(peaks_site,fimo_site_xw11,by=c("chr","start","end"), all.x = TRUE)
peaks_site <- merge(peaks_site,cut_site_xw11,by=c("chr","start","end"), all.x = TRUE)

n11 <- length(which(peaks_site$fimo & peaks_site$cut))
n10 <- length(which(peaks_site$fimo & is.na(peaks_site$cut)))
n01 <- length(which(is.na(peaks_site$fimo) & peaks_site$cut))
n00 <- length(which(is.na(peaks_site$fimo) & is.na(peaks_site$cut)))
fisher_result <- fisher.test(matrix(c(n11, n10, n01, n00), ncol = 2))
fisher_result
matrix(c(n11, n10, n01, n00), ncol = 2)

# bzip23
fimo_site <- read.table("/Users/lzz/Desktop/Rfiles/peak_unit/tf_motif.bed")

fimo_site_bzip23 <- fimo_site %>%
  filter(V7 == "ABF4") %>%
  distinct(V2,V3,.keep_all = TRUE) %>%
  setNames(c("chr","start","end")) %>%
  dplyr::select(1,2,3) %>%
  mutate(fimo = 1)

cut_site_bzip23 <- read.table("/Users/lzz/Desktop/Rfiles/cuttag/bzip23_RT.bed")
cut_site_bzip23 <- cut_site_bzip23 %>%
  distinct(V5,V6,.keep_all = TRUE) %>%
  dplyr::select(1,5,6) %>%
  setNames(c("chr","start","end")) %>%
  mutate(cut = 1)

peaks_site <- peaks_reads[,c(1,2,3)]
peaks_site <- merge(peaks_site,fimo_site_bzip23,by=c("chr","start","end"), all.x = TRUE)
peaks_site <- merge(peaks_site,cut_site_bzip23,by=c("chr","start","end"), all.x = TRUE)

n11 <- length(which(peaks_site$fimo & peaks_site$cut))
n10 <- length(which(peaks_site$fimo & is.na(peaks_site$cut)))
n01 <- length(which(is.na(peaks_site$fimo) & peaks_site$cut))
n00 <- length(which(is.na(peaks_site$fimo) & is.na(peaks_site$cut)))
fisher_result <- fisher.test(matrix(c(n11, n10, n01, n00), ncol = 2))
fisher_result
matrix(c(n11, n10, n01, n00), ncol = 2)
