setwd("~/Desktop/cuttag/chipqc/")
library(ChIPQC)

samplesheet <- read.csv("../samplesheet.csv")
chipObj <- ChIPQC(samplesheet, chromosomes = NULL) 
ChIPQCreport(chipObj)
