setwd("~/Desktop/cuttag/chipqc/")
library(ChIPQC)

samplesheet <- read.csv("../sampleSheet.csv")
chipObj <- ChIPQC(samplesheet, chromosomes = NULL) 
ChIPQCreport(chipObj)

