setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

gtsp <- pull.genoprob(mcross)

vita1a <- "1_3010272"

mp <- gtsp[, grep(vita1a, colnames(gtsp))]
gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
  if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
}))

longevity = as.numeric(pull.pheno(mcross)[, "longevity"])

ii1 <- which(gts == "AC")
ii2 <- which(gts == "AD")
ii3 <- which(gts == "BC")
ii4 <- which(gts == "BD")

col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")

nnCH <- length(ii1)
nnBH <- length(ii2)
nnCD <- length(ii3)
nnBD <- length(ii4)
#
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny")
pdf("KM-Vita1A-Combined.insetBottom.pdf")
par(cex=1.5)
par(cex.axis=1.2)
plot(c(900, 1200), c(50, 0), t = "n", main = "Vita1a -Combined - % KM Curve", ylab = "% of animals remaining", xlab = "Age", yaxt = "n")
mm <- c()
for(x in seq(15, 1500, 15)){
  iiD <- which(longevity >= x-15 & longevity < x)
  nCH <- length(which(iiD %in% ii1))
  nBH <- length(which(iiD %in% ii2))
  nCD <- length(which(iiD %in% ii3))
  nBD <- length(which(iiD %in% ii4))
  cat(x, "", nCH, nBH, nCD, nBD, "\n")
  mm <- rbind(mm, c(nCH, nBH, nCD, nBD))
  #points(x, nAC, col = 1)
  #points(x, nAD, col = 2)
  #points(x, nBC, col = 3)
  #points(x, nBD, col = 4)
}
colnames(mm) <- c("CH", "BH", "CD", "BD")
rownames(mm) <- seq(15, 1500, 15)
write.table(mm, "deathsIn15Dwindows.txt", sep = "\t", quote = FALSE)
axis(2, at = seq(0, 100, 2.5), seq(0, 100, 2.5), las=2)
points(seq(15, 1500, 15), 100 - (100 *cumsum(mm[,1]/nnCH)), t = "l", col = col.main[1])
points(seq(15, 1500, 15), 100 - (100 *cumsum(mm[,2]/nnBH)), t = "l", col = col.main[2])
points(seq(15, 1500, 15), 100 - (100 *cumsum(mm[,3]/nnCD)), t = "l", col = col.main[3])
points(seq(15, 1500, 15), 100 - (100 *cumsum(mm[,4]/nnBD)), t = "l", col = col.main[4])
dev.off()




#
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny")
#pdf("KM-Vita1A-DeathP15DWindow.pdf")
#par(cex=1.5)
#par(cex.axis=1.2)
plot(c(0, 1200), c(0.1, 50), t = "n", main = "Vita1a -Combined - % Deaths per Window", 
  ylab = "% of animals dead", xlab = "Age", yaxt = "n", log = "y")
mm <- c()

X <- seq(15, 1200, 15)

for(x in X){
  iiD <- which(longevity >= x-15 & longevity < x)
  nCH <- length(which(iiD %in% ii1))
  nBH <- length(which(iiD %in% ii2))
  nCD <- length(which(iiD %in% ii3))
  nBD <- length(which(iiD %in% ii4))
  cat(x, "", nCH, nBH, nCD, nBD, "\n")
  mm <- rbind(mm, c(nCH, nBH, nCD, nBD))
  #points(x, nAC, col = 1)
  #points(x, nAD, col = 2)
  #points(x, nBC, col = 3)
  #points(x, nBD, col = 4)
}
colnames(mm) <- c("CH", "BH", "CD", "BD")
rownames(mm) <- X
write.table(mm, "deathsIn15Dwindows.txt", sep = "\t", quote = FALSE)
Y1 <- 100 *mm[,1] / (nnCH - cumsum(mm[,1]))
#points(X, predict(loess(Y1 ~ X + 0)), t = "l", col = col.main[1])
points(X, Y1, t = "p", col = col.main[1])
abline(lm(log(Y1) ~ X), col = col.main[1])
Y2 <- 100 *mm[,2]/ (nnBH - cumsum(mm[,2]))
points(X, Y2, t = "p", col = col.main[2])
#points(X, predict(loess(Y2 ~ X + 0)), t = "l", col = col.main[2])
Y3 <- 100 *mm[,3]/ (nnCD- cumsum(mm[,3]))
points(X, Y3, t = "p", col = col.main[3])
#points(X, predict(loess(Y3 ~ X + 0)), t = "l", col = col.main[3])
Y4 <- 100 *mm[,4]/ (nnBD- cumsum(mm[,4]))
points(X, Y4, t = "p", col = col.main[4])
#points(X, predict(loess(Y4 ~ X + 0)), t = "l", col = col.main[4])
#dev.off()
axis(2, at = seq(0, 50, 2.5), seq(0, 50, 2.5),las=2)



nnCH <- length(ii1)
nnBH <- length(ii2)
nnCD <- length(ii3)
nnBD <- length(ii4)


pCH <- nnCH / 20 
pBH <- nnBH / 20 
pCD <- nnCD / 20 
pBD <- nnBD / 20 

round(seq(0, nnCH, pCH))
#setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny")
#pdf("TimeFor5percToDie.pdf")
#par(cex=1.5)
#par(cex.axis=1.2)
plot(c(1, 20), c(10, 500), t = "n", log = "y", main = "Time for 5% of the population to Die", 
     xaxt = "n", yaxt = "n", xlab = "5% mortality window", ylab = "Mortality window duration (d)")
points(diff(c(0, sort(longevity[ii1])[round(seq(0, nnCH, pCH))])), t = "l", col = col.main[1])
points(diff(c(0, sort(longevity[ii2])[round(seq(0, nnBH, pBH))])), t = "l", col = col.main[2])
points(diff(c(0, sort(longevity[ii3])[round(seq(0, nnCD, pCD))])), t = "l", col = col.main[3])
points(diff(c(0, sort(longevity[ii4])[round(seq(0, nnBD, pBD))])), t = "l", col = col.main[4])
axis(1, at = seq(1,20,1), seq(5, 100, 5))
axis(2, at = seq(0,500,10), seq(0, 500, 10), las=2)

#dev.off()
