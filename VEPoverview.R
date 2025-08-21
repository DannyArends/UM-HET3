#
# VEPoverview.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Print an overview of the number of MODERATE and HIGH impact SNPs per chromosome (VEP predition)
#

setwd("/home/rqdt9/Data/UM-HET3/wholeGenome")

hasM <- c()
hasH <- c()

for(chr in c(1:19, "X")){
  incon <- gzfile(paste0("chr",chr,".vep.gz"))
  mdata <- readLines(incon)
  mdata <- mdata[-c(1:80)]

  for(x in 1:length(mdata)){
    isM <- grepl("MODERATE", mdata[x])
    isH <- grepl("HIGH", mdata[x])
    if(isM) hasM <- c(hasM, strsplit(mdata[x], "\t")[[1]][4])
    if(isH) hasH <- c(hasH, strsplit(mdata[x], "\t")[[1]][4])
  }
  close(incon)
  cat("Done Chr", chr, "nModerate=", length(unique(hasM)), "nHigh=",length(unique(hasH)), "\n")
}

