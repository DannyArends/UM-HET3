setwd("C:/Users/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Users/rqdt9/OneDrive - Northumbria University - Production Azure AD/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 1)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

# Our Progressive Mapping Sequence
msequence <- seq(365, 1100, 15)
markers <- unique(unlist(lapply(strsplit(colnames(gtsp), ":"), "[",1)))

lods.cM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] <= x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  minAge <- c(minAge, min(cdata[, "longevity"]))

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
  cat("Done", x, "\n")
}
colnames(lods.cM) <- colnames(pull.geno(mcross))

png("ProgressiveMapping_UMHET3_DieBefore.png", width = 1400, height = 1050)

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "Purples")

layout(matrix(c(1,1,1,1,1,1,2,3,3,3,3,3,3,4,5,5,5,5,5,5,6), ncol=7, byrow=TRUE))

op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.cM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females and males", breaks = c(0, 2, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
#contour(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.cM), levels = c(0, 3.65, 4.25, 4.95, 5.95, 100), add = TRUE)

abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0("<=", msequence, " days")
axis(2, at = (1:length(yaxisD))[seq(1,length(yaxisD), 4)], yaxisD[seq(1,length(yaxisD), 4)], las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
box()
op <- par(mar = c(4,0,3,0))
plot(1:10, t = 'n', bty="n", xaxt ='n', yaxt ='n',xlab="")
legend("top", legend= c("<2.00", "2.00 - 3.65", "3.65 - 4.25", "4.25 - 4.95", "4.95 - 5.95", ">5.95"), fill = c("white", colz), bg="white")


### females

lods.fM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] <= x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  females <- which(cdata[, "sex"] == 0)
  cdata <- cdata[females,]
  gtsM <- gtsM[females, ]
  
  minAge <- c(minAge, min(cdata[, "longevity"]))

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
  cat("Done", x, "\n")
}
colnames(lods.fM) <- colnames(pull.geno(mcross))

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "PuRd")[-1]

op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.fM), 1:nrow(lods.fM), t(lods.fM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females", breaks = c(0, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0("<=", msequence, " days")
axis(2, at = (1:length(yaxisD))[seq(1,length(yaxisD), 4)], yaxisD[seq(1,length(yaxisD), 4)], las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
box()
op <- par(mar = c(4,0,3,0))
plot(1:10, t = 'n', bty="n", xaxt ='n', yaxt ='n',xlab="")
legend("top", legend= c("<3.65", "3.65 - 4.25", "4.25 - 4.95", "4.95 - 5.95", ">5.95"), fill = c("white", colz), bg="white")


### males
lods.mM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] <= x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  males <- which(cdata[, "sex"] == 1)
  cdata <- cdata[males,]
  gtsM <- gtsM[males, ]
  
  minAge <- c(minAge, min(cdata[, "longevity"]))

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
  cat("Done", x, "\n")
}
colnames(lods.mM) <- colnames(pull.geno(mcross))

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "Blues")

op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.mM), 1:nrow(lods.mM), t(lods.mM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - males", breaks = c(0,2, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0("<=", msequence, " days")
axis(2, at = (1:length(yaxisD))[seq(1,length(yaxisD), 4)], yaxisD[seq(1,length(yaxisD), 4)], las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
box()
op <- par(mar = c(4,0,3,0))
plot(1:10, t = 'n', bty="n", xaxt ='n', yaxt ='n',xlab="")
legend("top", legend= c("<3.65", "3.65 - 4.25", "4.25 - 4.95", "4.95 - 5.95", ">5.95"), fill = c("white", colz), bg="white")

dev.off()

rownames(lods.mM) <- paste0("> ", msequence)
rownames(lods.fM) <- paste0("> ", msequence)
rownames(lods.cM) <- paste0("> ", msequence)

write.table(round(lods.mM,2), "progressiveMapping_males_diebefore.txt", sep = "\t", quote=FALSE)
write.table(round(lods.fM,2), "progressiveMapping_females_diebefore.txt", sep = "\t", quote=FALSE)
write.table(round(lods.cM,2), "progressiveMapping_all_diebefore.txt", sep = "\t", quote=FALSE)

