# fisher test
library(dplyr)

# x39
peaks_reads <- read.table("~/Desktop/Rfiles/idr_peaks/peaks_reads_2_control.txt", header = T, sep = "\t")
fimo_site <- read.table("/Users/lzz/Desktop/Rfiles/peak_unit/motif_RT.bed", sep = "\t")
colnames(fimo_site)[8] <- "ID"
symbol <- read.table("/Users/lzz/Desktop/Rfiles/symbol.txt")
colnames(symbol) <- c("ID","TF")
fimo_site <- merge(fimo_site,symbol,by="ID")
fimo_site <- fimo_site[,c(2,3,4,6,7,8,14,1,9,13,5)]

fimo_site_x39 <- fimo_site %>%
  filter(TF == "ATHB-40") %>%
  distinct(V2,V3,.keep_all = TRUE) %>%
  dplyr::select(1,2,3) %>%
  setNames(c("chr","start","end")) %>%
  mutate(fimo = 1)

cut_site_x39 <- read.table("/Users/lzz/Desktop/Rfiles/cuttag/x39_RT.bed")
cut_site_x39 <- cut_site_x39 %>%
  distinct(V5,V6,.keep_all = TRUE) %>%
  dplyr::select(1,2,3) %>%
  setNames(c("chr","start","end")) %>%
  mutate(cut = 1)

peaks_site <- peaks_reads[,c(1,2,3,11)]
peaks_site <- merge(peaks_site,fimo_site_x39,by=c("chr","start","end"), all.x = TRUE)
peaks_site <- merge(peaks_site,cut_site_x39,by=c("chr","start","end"), all.x = TRUE)

n11_E <- length(which(peaks_site$fimo & peaks_site$cut & peaks_site$RT == "E"))
n10_E <- length(which(peaks_site$fimo & is.na(peaks_site$cut) & peaks_site$RT == "E"))
n01_E <- length(which(is.na(peaks_site$fimo) & peaks_site$cut & peaks_site$RT == "E"))
n00_E <- length(which(is.na(peaks_site$fimo) & is.na(peaks_site$cut) & peaks_site$RT == "E"))
E1 <- n11_E+n00_E ; E0 <- n10_E+n01_E

n11_M <- length(which(peaks_site$fimo & peaks_site$cut & peaks_site$RT == "M"))
n10_M <- length(which(peaks_site$fimo & is.na(peaks_site$cut) & peaks_site$RT == "M"))
n01_M <- length(which(is.na(peaks_site$fimo) & peaks_site$cut & peaks_site$RT == "M"))
n00_M <- length(which(is.na(peaks_site$fimo) & is.na(peaks_site$cut) & peaks_site$RT == "M"))
M1 <- n11_M+n00_M ; M0 <- n10_M+n01_M

n11_L <- length(which(peaks_site$fimo & peaks_site$cut & peaks_site$RT == "L"))
n10_L <- length(which(peaks_site$fimo & is.na(peaks_site$cut) & peaks_site$RT == "L"))
n01_L <- length(which(is.na(peaks_site$fimo) & peaks_site$cut & peaks_site$RT == "L"))
n00_L <- length(which(is.na(peaks_site$fimo) & is.na(peaks_site$cut) & peaks_site$RT == "L"))
L1 <- n11_L+n00_L ; L0 <- n10_L+n01_L

n11_o <- length(which(peaks_site$fimo & peaks_site$cut & !(peaks_site$RT %in% c("E","M","L"))))
n10_o <- length(which(peaks_site$fimo & is.na(peaks_site$cut) & !(peaks_site$RT %in% c("E","M","L"))))
n01_o <- length(which(is.na(peaks_site$fimo) & peaks_site$cut & !(peaks_site$RT %in% c("E","M","L"))))
n00_o <- length(which(is.na(peaks_site$fimo) & is.na(peaks_site$cut) & !(peaks_site$RT %in% c("E","M","L"))))
O1 <- n11_o+n00_o ; O0 <- n10_o+n01_o

fisher_result <- chisq.test(matrix(c(E1, E0, M1, M0, L1, L0, O1, O0), ncol = 4))
fisher_result
matrix(c(n11, n10, n01, n00), ncol = 2)


# x49
peaks_reads <- read.table("~/Desktop/Rfiles/idr_peaks/peaks_reads_2_control.txt", header = T, sep = "\t")
fimo_site <- read.table("/Users/lzz/Desktop/Rfiles/peak_unit/motif_RT.bed", sep = "\t")
colnames(fimo_site)[8] <- "ID"
symbol <- read.table("/Users/lzz/Desktop/Rfiles/symbol.txt")
colnames(symbol) <- c("ID","TF")
fimo_site <- merge(fimo_site,symbol,by="ID")
fimo_site <- fimo_site[,c(2,3,4,6,7,8,14,1,9,13)]

fimo_site_x49 <- fimo_site %>%
  filter(TF == "HAT5") %>%
  distinct(V2,V3,.keep_all = TRUE) %>%
  setNames(c("chr","start","end")) %>%
  dplyr::select(1,2,3) %>%
  mutate(fimo = 1)

cut_site_x49 <- read.table("/Users/lzz/Desktop/Rfiles/cuttag/x49_RT.bed")
cut_site_x49 <- cut_site_x49 %>%
  distinct(V5,V6,.keep_all = TRUE) %>%
  dplyr::select(1,2,3) %>%
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
peaks_reads <- read.table("~/Desktop/Rfiles/idr_peaks/peaks_reads_2_control.txt", header = T, sep = "\t")
fimo_site <- read.table("/Users/lzz/Desktop/Rfiles/peak_unit/motif_RT.bed", sep = "\t")
colnames(fimo_site)[8] <- "ID"
symbol <- read.table("/Users/lzz/Desktop/Rfiles/symbol.txt")
colnames(symbol) <- c("ID","TF")
fimo_site <- merge(fimo_site,symbol,by="ID")
fimo_site <- fimo_site[,c(2,3,4,6,7,8,14,1,9,13)]

fimo_site_xw11 <- fimo_site %>%
  filter(TF == "WOX11") %>%
  distinct(V2,V3,.keep_all = TRUE) %>%
  setNames(c("chr","start","end")) %>%
  dplyr::select(1,2,3) %>%
  mutate(fimo = 1)

cut_site_xw11 <- read.table("/Users/lzz/Desktop/Rfiles/cuttag/xw11_RT.bed")
cut_site_xw11 <- cut_site_xw11 %>%
  distinct(V5,V6,.keep_all = TRUE) %>%
  dplyr::select(1,2,3) %>%
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

