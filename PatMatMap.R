library(svglite)

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
mapAll <- scanone(mcross)[,1:2]

gMat <- read.table(file="ITP_6480x486_Maternal_Sep21.txt", sep = "\t",na.strings=c("", "NA", "x", "Not available"))
gMap <- gMat[,1:2]
gMat <- gMat[,-c(1:2)]
gPat <- read.table(file="ITP_6480x396_Paternal_Sep21.txt", sep = "\t",na.strings=c("", "NA", "x", "Not available"))
gPap <- gPat[,1:2]
gPat <- gPat[,-c(1:2)]


setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/01_Paper_and_Main_Arends_Working_Files_2023/Supplemental files")
svglite(paste0("geneticmap.svg"), width = 24, height = 12)


plot(x = c(1,20), y = c(-2, 195), t = 'n', xlab = "Chromosome", ylab = "Location (Mb)", xaxt='n', yaxt='n', las=2, yaxs='i', main = "UM-HET3 genetic marker map")
cnt <- 1
for(chr in c(1:19,"X")) {
  aa <- rownames(mapAll)[which(mapAll[,1] == chr)]
  colz <- rep("black", length(aa))
  colz[which(aa %in% rownames(gMap))] <- "hotpink"
  colz[which(aa %in% rownames(gPap))] <- "blue"
  for(x in 1:length(aa)){
    lines(c(cnt - 0.2, cnt + 0.2), c(mapAll[aa[x],2],mapAll[aa[x],2]), col=colz[x])
  }
  rect(cnt - 0.2, 0, cnt + 0.2, max(mapAll[aa,2]))
  cnt <- cnt + 1
}
axis(1, at = 1:20, c(1:19, "X"))
axis(2, at = seq(0,200,25), seq(0,200,25),las=2)
legend("topright", c("Maternal", "Paternal", "X chromosome"), lwd=1, col = c("hotpink", "blue", "black"))

dev.off()



aa <- cor(pull.geno(fill.geno(mcross)), use = "pair")

library(RColorBrewer)
colz <- colorRampPalette(c("white", "blue", "green"))(50)
image(1:nrow(aa), 1:ncol(aa), abs(aa), breaks = seq(0,1,0.02), col=colz, main = "Correlation between markers", xlab="Chromosome", xaxt='n', ylab="Chromosome", yaxt='n')

abline(h = cumsum(table(mapAll[,1])))
abline(v = cumsum(table(mapAll[,1])))
axis(1, at = (c(0, cumsum(table(mapAll[,1]))) +  0.5 * diff(c(0, cumsum(table(mapAll[,1])))))[-21], c(1:19,"X"),cex.axis=0.8)
axis(2, at = (c(0, cumsum(table(mapAll[,1]))) +  0.5 * diff(c(0, cumsum(table(mapAll[,1])))))[-21], c(1:19,"X"), las=2,cex.axis=0.8)


#

annot <- read.table("phenotypes/HET3-ITP_traitsGN2_Sep2021.csv", sep=",", skip=11, header=TRUE,na.strings=c("x", "Not available"))
annot <- annot[-c(1:20),]
annot <- annot[-grep("_SE", annot[,1]),]
rownames(annot) <- annot[,1]
annot <- annot[,-1]

# Animals marked as REMOVED, but genotyped
removed <- c("UT04802", "UT04803", "UT04501", "UT04502", "UT04539", "UT04540", 
             "UT04541", "UT04542", "UT04725", "UT04726", "UT04727", "UT04728",
             "UM39151", "UM40061", "UM40467", "UM41174", "UM41419", "UM41426", "UM41993", "UM42007", "UM43218")

annot <- cbind(annot, sex = 1 - as.numeric(annot[,"Sex"]))
annot <- annot[-which(annot[, "Tx_Group"] == "2.0"),]
annot <- annot[-which(rownames(annot) %in% removed),]


fdata <- data.frame(longevity = annot[, "Longevity_HET3_ITP"], 
                    bw6 = annot[, "ITP_Weight6m"], 
                    bw12 = annot[, "ITP_Weight12m"], 
                    bw18 = annot[, "Body_18m"], 
                    bw24 = annot[, "BodyWeight_HET3_ITP_24m"], 
                    sex = as.factor(annot[, "sex"]), 
                    site = as.factor(annot[, "Site"]),
                    cohort = as.factor(annot[, "Year"]), 
                    treatment = as.factor(annot[, "Tx_Group"]))
rownames(fdata) <- rownames(annot)
idx <- rownames(fdata)[which(fdata[, "longevity"] > 365)]
idx <- idx[which(idx %in% colnames(gMat) & idx %in% colnames(gPat))]

fdata <- fdata[idx,]

gMat <- gMat[,idx]
gPat <- gPat[,idx]

msequence <- seq(0, 0.8, 0.02)

## Maternal Map

## Combined longevity > 356 days
lods.cM <- c()
minAge <- c()
for(x in msequence){
  lods.c <- c()
  idxs <- sort(fdata[, "longevity"], index.return=TRUE)$ix
  idxs <- idxs[(1+ (x * length(idxs))): length(idxs)]
  cdata <- fdata[idxs,]
  minAge <- c(minAge, min(cdata[, "longevity"]))
  gMatS <- gMat[, idxs]
  cat("dim:", dim(cdata)[1], "\n")
  for(marker in rownames(gMatS)){
    mp <- as.character(as.numeric(gMatS[marker,]))
    tryCatch(
    {
      lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp, data = cdata)
      lods.c <<- c(lods.c, -log10(anova(lm.alt)["mp", "Pr(>F)"]))
    },error = function(cond) {
      lods.c <<- c(lods.c, NA)
    })
  }
  names(lods.c) <- rownames(gMat)
  lods.cM <- rbind(lods.cM, lods.c)
}
rownames(lods.cM) <- msequence

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(4, "BuGn")
op <- par(mar = c(4,10,3,1))
image(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.cM), xlab="Maternal map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females and males", breaks = c(0, 3, 4, 6, 8, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(gMap[,"Chr"]))) != 0))
yaxisD <- paste0(c("All (>", paste("top", 100*(1 - msequence[-1]), "% (>")), minAge, " days)")
axis(2, at = 1:length(yaxisD), yaxisD, las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(gMap[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(gMap)))/2), c(1:19))
legend("topright", legend= c("0-3", "3-4", "4-6", "6-8", "8+"), fill = c("white", colz), bg="white")
box()

## Paternal Map

## Combined longevity > 356 days
lods.cP <- c()
minAge <- c()
for(x in msequence){
  lods.c <- c()
  idxs <- sort(fdata[, "longevity"], index.return=TRUE)$ix
  idxs <- idxs[(1+ (x * length(idxs))): length(idxs)]
  cdata <- fdata[idxs,]
  minAge <- c(minAge, min(cdata[, "longevity"]))
  gPatS <- gPat[, idxs]
  cat("dim:", dim(cdata)[1], "\n")
  for(marker in rownames(gPatS)){
    mp <- as.character(as.numeric(gPatS[marker,]))
    tryCatch(
    {
      lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp, data = cdata)
      lods.c <<- c(lods.c, -log10(anova(lm.alt)["mp", "Pr(>F)"]))
    },error = function(cond) {
      lods.c <<- c(lods.c, NA)
    })
  }
  names(lods.c) <- rownames(gPat)
  lods.cP <- rbind(lods.cP, lods.c)
}
rownames(lods.cP) <- msequence


# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(4, "BuGn")
op <- par(mar = c(4,10,3,1))
image(1:ncol(lods.cP), 1:nrow(lods.cP), t(lods.cP), xlab="Maternal map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females and males", breaks = c(0, 3, 4, 6, 8, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(gMap[,"Chr"]))) != 0))
yaxisD <- paste0(c("All (>", paste("top", 100*(1 - msequence[-1]), "% (>")), minAge, " days)")
axis(2, at = 1:length(yaxisD), yaxisD, las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(gMap[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(gMap)))/2), c(1:19))
legend("topright", legend= c("0-3", "3-4", "4-6", "6-8", "8+"), fill = c("white", colz), bg="white")
box()

apply(lods.cM,1,function(x){max(x,na.rm=TRUE)})
apply(lods.cP,1,function(x){max(x,na.rm=TRUE)})



which(apply(lods.cM,1,function(x){max(x,na.rm=TRUE)}) > 4.25)
