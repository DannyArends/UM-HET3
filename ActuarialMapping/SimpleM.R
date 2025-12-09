#
# simpleM.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Use the simpleM method to iteratively determine how many independant tests are performed when performing a single genotype scan
# see: https://github.com/LTibbs/SimpleM
#

setwd("C:/Users/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Users/rqdt9/OneDrive - Northumbria University - Production Azure AD/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))

days <- as.numeric(pull.pheno(mcross)[, "longevity"])
idx <- which(days >= 365)
mgeno <- pull.geno(fill.geno(mcross))
numgeno <- t(mgeno[idx,])

Meff_PCA <- function(eigenValues, PCA_cutoff = 0.995){
  totalEigenValues <- sum(eigenValues)
  myCut <- PCA_cutoff*totalEigenValues
  num_Eigens <- length(eigenValues)
  myEigenSum <- 0
  index_Eigen <- 0
  
  for(i in 1:num_Eigens){
    if (myEigenSum <= myCut) {
      myEigenSum <- myEigenSum + eigenValues[i]
      index_Eigen <- i
    } else {
      break
    }
  }  
  return(index_Eigen)
}

inferCutoff <- function(dt_My, PCA_cutoff = 0.995){ # infer the cutoff => Meff
  CLD <- cor(dt_My, use = "pair")
  eigen_My <- eigen(CLD)
  eigenValues_dt <- abs(eigen_My$values)
  Meff_PCA_gao <- Meff_PCA(eigenValues_dt, PCA_cutoff)
  return(Meff_PCA_gao)
}

numLoci <- nrow(numgeno)
simpleMeffs <- NULL
rangeTested <- seq(10, 500, 10)
for(x in rangeTested){
  simpleMeff <- NULL
  fixLength <- x 
  i <- 1
  myStart <- 1
  myStop <- 1
  while (myStop < numLoci) {
    myDiff <- numLoci - myStop 
    if (myDiff <= fixLength) break
    
    myStop <- myStart + i*fixLength - 1
    snpInBlk <- t(numgeno[myStart:myStop, ])
    MeffBlk <- inferCutoff(snpInBlk)
    simpleMeff <- c(simpleMeff, MeffBlk)
    myStart <- myStop+1
  }
  snpInBlk <- t(numgeno[myStart:numLoci, ])
  if (nrow(snpInBlk) > 1) {
    MeffBlk <- inferCutoff(snpInBlk)
  } else {
    MeffBlk <- 1
  }
  simpleMeff <- c(simpleMeff, MeffBlk)

  cat("Total number of SNPs is: ", numLoci, "\n")
  cat("Inferred Meff is: ", sum(simpleMeff), "\n")
  simpleMeffs <- c(simpleMeffs, sum(simpleMeff))
}

c(fixLength = rangeTested[which.min(simpleMeffs)], nTests = min(simpleMeffs))
plot(rangeTested, simpleMeffs, xlab = "fixLength", ylab="Meff")
