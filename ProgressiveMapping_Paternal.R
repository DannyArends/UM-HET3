#
# ProgressiveMapping_Maternal.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Map longevity QTLs by progressively raising the threshold for inclusion (T-age) by 15 days (starting at 20 days)
# We perform a full scan at each T-age for both Combined (F+M)), Females, and Males
# We don't use 4 genotype probabilities at each marker, but use paternal haplotypes (C and D)
#

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

# Our Progressive Mapping Sequence
msequence <- seq(365, 1100, 15)

lods.cM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  # Subset data based on the current threshold
  idx <- which(cdata[, "longevity"] >= x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  lods.c <- c()
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
      if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
    }))

    gts[gts == "AC"] <- "C"; gts[gts == "AD"] <- "D"
    gts[gts == "BC"] <- "C"; gts[gts == "BD"] <- "D"

    iix <- which(!is.na(gts))
    gts <- gts[iix]
    sdata <- cdata[iix, ]

    tryCatch(
      {
      lm.null <- lm(longevity ~ sex + site + cohort + treatment + 0, data = sdata)
      lm.alt <- lm(longevity ~ sex + site + cohort + treatment + gts + 0, data = sdata)
      n <- sum(!is.na(lm.alt$resid))
      lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
      lods.c <- c(lods.c, lod)
      },error=function(cond) {
        lods.c <<- c(lods.c, NA)
      }
    )
  }
  cat("Done",x, "\n")
  lods.cM <- rbind(lods.cM, lods.c)
}
colnames(lods.cM) <- colnames(pull.geno(mcross))
rownames(lods.cM) <- paste0("> ", msequence)

write.table(round(lods.cM,2), "progressiveMapping_pat_all.txt", sep = "\t", quote=FALSE)

### Females

lods.fM <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  # Subset data based on the current threshold
  idx <- which(cdata[, "longevity"] >= x & cdata[, "sex"] == 0)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  lods.c <- c()

  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
      if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
    }))

    gts[gts == "AC"] <- "C"; gts[gts == "AD"] <- "D"
    gts[gts == "BC"] <- "C"; gts[gts == "BD"] <- "D"

    iix <- which(!is.na(gts))
    gts <- gts[iix]
    sdata <- cdata[iix, ]
    tryCatch(
      {
      lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = sdata)
      lm.alt <- lm(longevity ~ site + cohort + treatment + gts + 0, data = sdata)
      n <- sum(!is.na(lm.alt$resid))
      lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
      lods.c <- c(lods.c, lod)
      },error=function(cond) {
        lods.c <<- c(lods.c, NA)
      }
    )
  }
  cat("Done",x, "\n")
  lods.fM <- rbind(lods.fM, lods.c)
}

colnames(lods.fM) <- colnames(pull.geno(mcross))
rownames(lods.fM) <- paste0("> ", msequence)

write.table(round(lods.fM,2), "progressiveMapping_pat_females.txt", sep = "\t", quote=FALSE)

### Males

lods.mM <- c()
for(x in seq(365, 1050, 15)){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  # Subset data based on the current threshold
  idx <- which(cdata[, "longevity"] >= x & cdata[, "sex"] == 1)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  lods.c <- c()

  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
      if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
    }))

    gts[gts == "AC"] <- "C"; gts[gts == "AD"] <- "D"
    gts[gts == "BC"] <- "C"; gts[gts == "BD"] <- "D"

    iix <- which(!is.na(gts))
    gts <- gts[iix]
    sdata <- cdata[iix, ]
    tryCatch(
      {
      lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = sdata)
      lm.alt <- lm(longevity ~ site + cohort + treatment + gts + 0, data = sdata)
      n <- sum(!is.na(lm.alt$resid))
      lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
      lods.c <- c(lods.c, lod)
      },error=function(cond) {
        lods.c <<- c(lods.c, NA)
      }
    )
  }
  cat("Done",x, "\n")
  lods.mM <- rbind(lods.mM, lods.c)
}

colnames(lods.mM) <- colnames(pull.geno(mcross))
rownames(lods.mM) <- paste0("> ", seq(365, 1050, 15))

write.table(round(lods.mM,2), "progressiveMapping_pat_males.txt", sep = "\t", quote=FALSE)


