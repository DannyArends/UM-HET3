#
# Figure1.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Figures for the main paper text (Main text Figure 1, but some figure might have moved into other figures to make them fit)
# Figures: 
# - Kaplan–Meier curve (male & female) (1b)
# - Our Effect curve (Difference in Expectancy (days)) (1b)
# - Histograms of survivors both horizontal as well as vertical (Fig. 1a)
#

library(svglite)
library(pspline)
library(qtl)

source("ActuarialMapping/adjustXprobs.R")

# Read cross object
mcross <- read.cross(format="csvr", file="DataSet/um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

phe <- pull.pheno(mcross)[, "longevity"]
sex <- pull.pheno(mcross)[, "sex"]

males <- phe[sex == 1 & phe >= 20]
females <- phe[sex == 0 & phe >= 20]

pdf("DataSet/output/Figure_1_KM.pdf", width = 16, height = 12)
par(cex=2)
par(cex.axis=1.2)
plot(c(20, 1456), c(0, 100), t = "n", ylab = "% survival", xlab = "days", yaxt="n", main = "KM Curve", xaxt = "n")
for(x in 20:1456){
  n.n <- length(which(males >= x))
  n.r <- length(which(females >= x))
  points(x, (n.n/length(males))* 100, pch=18, col = "blue")
  points(x, (n.r/length(females))* 100, pch=18, col = "pink")
}
axis(2, at = c(0,25,50,75,100), c(0,25,50,75,100), las = 2)
axis(1, at = seq(0, 1500, 100), seq(0, 1500, 100), las = 2)
legend("topright", c("Males", "Females"), col = c("blue", "hotpink"), pch=18)
dev.off()

pdf("DataSet/output/Figure_1_EC.pdf", width = 12, height = 12)
par(cex=2)
par(cex.axis=1.2)

plot(c(42, 1100), c(-50, 50), t = "n", ylab = "Difference in Expectancy (days)", xlab = "Lifespan Cut-off Age (days)", yaxt="n", main = "Effect Curve")
msequence <- seq(42, 1100, 15)
std <- function(x) sd(x)/sqrt(length(x))
mm <- c()
for(d in msequence){
  mAll <- mean(c(males[males >= d], females[females >= d]))
  mM <- mean(c(males[males >= d])) - mAll
  sM <- std(males[males >= d])
  mF <- mean(c(females[females >= d])) - mAll
  sF <- std(females[females >= d])
  mm <- rbind(mm, c(mM, sM, mF, sF))
}

# Combined
points(msequence, mm[,1], t = 'l', col = "blue", lwd=2)
polygon(c(msequence, rev(msequence)), c(mm[,1] + mm[,2], rev(mm[,1] - mm[,2])), col = rgb(0,0,1,0.5), border = NA)

points(msequence, mm[,3], t = 'l', col = "hotpink", lwd=2)
polygon(c(msequence, rev(msequence)), c(mm[,3] + mm[,4], rev(mm[,3] - mm[,4])), col = rgb(1, 0.41, 0.7,0.5), border = NA)

legend("topright", c("Males", "Females"), col = c("blue", "hotpink"), pch=18)
axis(2, at = seq(-100, 100, 10), seq(-100, 100, 10), las=2)
dev.off()

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
pdf("DataSet/output/Histogram_Horizontal_M+F.pdf")
plot(x= c(-125, 125), y = c(0, 1500), main = "UM-HET3 longevity by sex",
     ylab = "Longevity (days)", xlab = "Frequency", t = "n", yaxs= "i", xaxs= "i", yaxt="n")
for(x in 1:length(m$counts)){
  rect(-m$counts[x], seq(0, 1500, 15)[x], 0, seq(0, 1500, 15)[x + 1], col = "cornflowerblue")
  rect(0, seq(0, 1500, 15)[x], f$counts[x], seq(0, 1500, 15)[x + 1], col = "hotpink")
}
axis(2, at = seq(0, 1500, 250), seq(0, 1500, 250), las=2)
dev.off()

pdf("DataSet/output/Histogram_Horizontal_M+F_remaining.pdf")
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

