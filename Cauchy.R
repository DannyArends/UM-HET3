cauchycombination <- function(p) {
   qc <- qcauchy(p, lower.tail=FALSE)
   pc <- pcauchy(mean(qc), lower.tail=FALSE)
   pc
}

setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

lods <- read.table("progressiveMapping_all.txt", sep = "\t")


cc <- c()
for(x in 1:ncol(lods)){
  cc <- c(cc, cauchycombination(10 ^ -as.numeric(lods[,x])))
}
