
library(svglite)
library(pspline)


setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)


phe <- pull.pheno(mcross)[, "longevity"]
sex <- pull.pheno(mcross)[, "sex"]

males <- phe[sex == 1 & phe >= 365]
females <- phe[sex == 0 & phe >= 365]

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
svglite("Figure_1_KM.svg", width = 12, height = 12)
par(cex=2)
par(cex.axis=1.2)
plot(c(365, 1456), c(0, 100), t = "n", ylab = "% survival", xlab = "days", yaxt="n", main = "KM Curve")
for(x in 365:1456){
  n.n <- length(which(males >= x))
  n.r <- length(which(females >= x))
  points(x, (n.n/length(males))* 100, pch=18, col = "blue")
  points(x, (n.r/length(females))* 100, pch=18, col = "pink")
}
axis(2, at = c(0,25,50,75,100), c(0,25,50,75,100), las=2)
legend("topright", c("Males", "Females"), col = c("blue", "hotpink"), pch=18)
dev.off()

svglite("Figure_1_EC.svg", width = 12, height = 12)
par(cex=2)
par(cex.axis=1.2)

plot(c(365, 1100), c(-35, 35), t = "n", ylab = "Difference in Expectancy (days)", xlab = "Lifespan Cut-off Age (days)", yaxt="n", main = "Effect Curve")
msequence <- seq(365, 1100, 15)
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
axis(2, at = seq(-100, 100, 5), seq(-100, 100, 5), las=2)
dev.off()




