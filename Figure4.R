setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
pheno <- pull.pheno(mcross)


# Which rows contains males (isM) and which rows contain females(isF)
isM <- which(pheno[ , "sex"] == 1)
isF <- which(pheno[ , "sex"] == 0)

plot.new()
m <- hist(pheno[isM, "longevity"], add = TRUE, 
          col = rgb(0, 0, 1, 0.3), breaks = seq(20, 1500, 15))
f <- hist(pheno[isF, "longevity"], add = TRUE, 
          col = rgb(1, 0.3, 0.7, 0.3), breaks = seq(20, 1500, 15))

# Fancy Schmacy histogram (upto Sachini's standards)
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_Tables_Files/06_RandomFiguresForRob")
pdf("Histogram_Horizontal_M+F.pdf")
plot(x= c(-125, 125), y = c(0, 1500), main = "UM-HET3 longevity by sex",
     ylab = "Longevity (days)", xlab = "Frequency", t = "n", yaxs= "i", xaxs= "i", yaxt="n")
for(x in 1:length(m$counts)){
  rect(-m$counts[x], seq(0, 1500, 15)[x], 0, seq(0, 1500, 15)[x + 1], col = "cornflowerblue")
  rect(0, seq(0, 1500, 15)[x], f$counts[x], seq(0, 1500, 15)[x + 1], col = "hotpink")
}
axis(2, at = seq(0, 1500, 250), seq(0, 1500, 250), las=2)
dev.off()

pdf("Histogram_Horizontal_M+F_remaining.pdf")
plot(x= c(-125, 125), y = c(0, 1300), main = "",
     ylab = "Longevity (days)", xlab = "Frequency", t = "n", yaxs= "i", xaxs= "i", yaxt="n", xaxt="n")
for(x in 1:length(m$counts)){
  mcnt <- -100 * (sum(m$counts[x:length(m$counts)]) / sum(m$counts))
  fcnt <- 100 * (sum(f$counts[x:length(f$counts)]) / sum(f$counts))
  rect(mcnt, seq(0, 1500, 15)[x], 0, seq(0, 1500, 15)[x + 1], col = "cornflowerblue", border = "white")

  rect(0, seq(0, 1500, 15)[x], fcnt, seq(0, 1500, 15)[x + 1], col = "deeppink3", border = "white")
  rect(0, seq(0, 1500, 15)[x], -mcnt, seq(0, 1500, 15)[x + 1], col = "hotpink", border = "white")
}
axis(1, at = seq(-100, 100, 10), abs(seq(-100, 100, 10)), las=2)
axis(2, at = seq(0, 1500, 100), seq(0, 1500, 100), las=2)
axis(3, at = seq(-100, 0, 10), round(sum(m$counts) * (abs(seq(-100, 0, 10)) / 100), 0), las=2)
axis(3, at = seq(0, 100, 10), round(sum(f$counts) * (abs(seq(0, 100, 10)) / 100), 0), las=2)
dev.off()

