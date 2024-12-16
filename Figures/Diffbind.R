setwd("~/Desktop/repli-ATAC-seq/cuttag/chipqc/")
library(DiffBind)
samplesheet <- read.csv("../sampleSheet1.csv")
dbObj <- dba(sampleSheet = samplesheet)
dbObj_consensus <- dba.peakset(dbObj,
                                    consensus=c(DBA_FACTOR),
                                    minOverlap=2)
dbObj_consensus <- dba(dbObj_consensus,
                           mask=dbObj_consensus$masks$Consensus,
                           minOverlap=1)
dbObj_consensus
dba.plotVenn(dbObj, dbObj$masks$All)

########output the bind site########
consensus_peaks <- dba.peakset(dbObj_consensus, bRetrieve=TRUE)
# write.csv(consensus_peaks,"~/Desktop/cuttag/DiffBind/BindSite.csv", quote = F, row.names = F)

########Get the bind sites of each TF#######
library(dplyr)
bindsite <- read.csv("~/Desktop/cuttag/DiffBind/BindSite.csv")
bindsite$midpoint <- round((bindsite$start + bindsite$end) / 2, 0)
bindsite$midpoint_1 <- bindsite$midpoint + 1
bs_x39 <- bindsite %>% 
  filter(x39 > 0) %>% 
  select(seqnames,midpoint,midpoint_1)
bs_x49 <- bindsite %>% 
  filter(x49 > 0) %>% 
  select(seqnames,midpoint,midpoint_1)
bs_xw11 <- bindsite %>% 
  filter(xw11 > 0) %>% 
  select(seqnames,midpoint,midpoint_1)
write.csv(bs_x39,"~/Desktop/cuttag/DiffBind/BS_x39.csv", quote = F, row.names = F, col.names = F)
write.csv(bs_x49,"~/Desktop/cuttag/DiffBind/BS_x49.csv", quote = F, row.names = F, col.names = NA)
write.csv(bs_xw11,"~/Desktop/cuttag/DiffBind/BS_xw11.csv", quote = F, row.names = F, col.names = NA)
