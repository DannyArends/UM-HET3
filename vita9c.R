setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

vita9c <- "9_104091597"

plot(c(365, 1100), c(-40, 40), t = 'n')

msequence <- seq(365, 1100, 15)

mm <- pull.geno(fill.geno(mcross))
phe <- pull.pheno(mcross)[, "longevity"]
marker <- mm[,vita9c]

for(x in msequence){
  plot(which(phe > x) ~ marker)
}


