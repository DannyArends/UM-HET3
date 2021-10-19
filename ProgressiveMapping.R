setwd("D:/Ddrive/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Github/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

msequence <- seq(0, 0.8, 0.02)

lods.cM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] > 365)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  idxs <- sort(cdata[, "longevity"], index.return=TRUE)$ix
  idxs <- idxs[(1+ (x * length(idxs))): length(idxs)]
  cdata <- cdata[idxs,]
  minAge <- c(minAge, min(cdata[, "longevity"]))
  gtsM <- gtsM[idxs,]

  lods.c <- c()
  lm.null <- lm(longevity ~ sex + site + cohort + treatment + 0, data = cdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.cM <- rbind(lods.cM, lods.c)
}
colnames(lods.cM) <- colnames(pull.geno(mcross))

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(4, "BuGn")

op <- par(mar = c(4,10,3,1))
image(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.cM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females and males", breaks = c(0,3,4,6,8,100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(c("All (>", paste("top", 100*(1 - msequence[-1]), "% (>")), minAge, " days)")
axis(2, at = 1:length(yaxisD), yaxisD, las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
box()
legend("topright", legend= c("0-3", "3-4", "4-6", "6-8", "8+"), fill = c("white", colz), bg="white")


### females

lods.fM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] > 365)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  females <- which(cdata[, "sex"] == 0)
  cdata <- cdata[females,]
  gtsM <- gtsM[females, ]
  
  idxs <- sort(cdata[, "longevity"], index.return=TRUE)$ix
  idxs <- idxs[(1+ (x * length(idxs))): length(idxs)]
  cdata <- cdata[idxs,]
  minAge <- c(minAge, min(cdata[, "longevity"]))
  gtsM <- gtsM[idxs,]

  lods.c <- c()
  lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = cdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(longevity ~ site + cohort + treatment + mp + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.fM <- rbind(lods.fM, lods.c)
}
colnames(lods.fM) <- colnames(pull.geno(mcross))

op <- par(mar = c(4,10,3,1))
image(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.fM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females", breaks = c(0,3,4,6,8,100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(c("All (>", paste("top", 100*(1 - msequence[-1]), "% (>")), minAge, " days)")
axis(2, at = 1:length(yaxisD), yaxisD, las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
box()
legend("topright", legend= c("0-3", "3-4", "4-6", "6-8", "8+"), fill = c("white", colz), bg="white")


### males
lods.mM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] > 365)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  males <- which(cdata[, "sex"] == 1)
  cdata <- cdata[males,]
  gtsM <- gtsM[males, ]
  
  idxs <- sort(cdata[, "longevity"], index.return=TRUE)$ix
  idxs <- idxs[(1+ (x * length(idxs))): length(idxs)]
  cdata <- cdata[idxs,]
  minAge <- c(minAge, min(cdata[, "longevity"]))
  gtsM <- gtsM[idxs,]

  lods.c <- c()
  lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = cdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(longevity ~ site + cohort + treatment + mp + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.mM <- rbind(lods.mM, lods.c)
}
colnames(lods.mM) <- colnames(pull.geno(mcross))

op <- par(mar = c(4,10,3,1))
image(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.mM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - males", breaks = c(0,3,4,6,8,100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(c("All (>", paste("top", 100*(1 - msequence[-1]), "% (>")), minAge, " days)")
axis(2, at = 1:length(yaxisD), yaxisD, las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
box()
legend("topright", legend= c("0-3", "3-4", "4-6", "6-8", "8+"), fill = c("white", colz), bg="white")


