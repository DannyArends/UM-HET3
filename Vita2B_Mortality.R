# TODO: Follow up with Arthus for the supplemental files
# TODO: Loess smoother
# TODO: Output a folder of these for all Vita / haplotype

marker <- "2_112255823"

library(svglite)

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)

gtsp <- pull.genoprob(mcross)

mp <- gtsp[, grep(marker, colnames(gtsp))]
gts <- unlist(lapply(lapply(lapply(apply(mp,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
  if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
}))

areC <- grep("C", gts)
areD <- grep("D", gts)

gtsH <- rep(NA, length(gts))
gtsH[areC] <- "C"
gtsH[areD] <- "D"

phe <- cbind(lifespan = pull.pheno(mcross)[, "longevity"], sex = pull.pheno(mcross)[, "sex"], gtsH)
phe <- phe[which(phe[, "sex"] == "1" & !is.na(phe[, "gtsH"])),]

C <- as.numeric(phe[phe[,3] == "C" & as.numeric(phe[,1]) >= 20, 1])
D <- as.numeric(phe[phe[,3] == "D" & as.numeric(phe[,1]) >= 20, 1])

hist(C, col = rgb(1,0,0,0.5), breaks = seq(0, 1500, 15))
hist(D, add = TRUE, col = rgb(0,1,0,0.5), breaks = seq(0, 1500, 15))


plot(c(200, 1100), c(-10, 10), t = "n", ylab = "% life gained at tAge", xlab = "Lifespan Cut-off Age (days)", yaxt="n", xaxs="i", main = "Effect Curve")
msequence <- seq(20, 1100, 15)
std <- function(x) sd(x)/sqrt(length(x))
mm <- c()
for(d in msequence){
  mAll <- mean(c(C, D))
  mM <- mean(c(C[C >= d])) - mAll
  sM <- mean(D[D >= d]) - mAll
  mm <- rbind(mm, c(mM, sM))
}
abline(h = seq(0, 15, 1), lty=2, col = "gray")
points(msequence, 100*((mm[,1]-mm[,2]) / ((mean(c(C,D)) * 2) - msequence)), t = 'l', col = "pink", lwd=2)
points(msequence, 100*((mm[,2]-mm[,1]) / ((mean(c(C,D)) * 2) - msequence)), t = 'l', col = "brown", lwd=2)
#polygon(c(msequence, rev(msequence)), c(mm[,1] + mm[,2], rev(mm[,1] - mm[,2])), col = rgb(0,0,1,0.5), border = NA)

#points(msequence, mm[,3], t = 'l', col = rgb(.4,0.2,1,1), lwd=2)
#polygon(c(msequence, rev(msequence)), c(mm[,3] + mm[,4], rev(mm[,3] - mm[,4])), col = rgb(.4,0.2,1,0.5), border = NA)


legend("topright", c("C relative to D"), col = c(rgb(0,0,1,0.5), rgb(1, 0.41, 0.7,0.5)), pch=18)
axis(2, at = seq(-100, 170, 1), seq(-100, 170, 1), las=2)
#dev.off()


setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny")
pdf("Vita2B_Deaths_15day_WholePopulation.pdf", width = 16, height = 12)

par(cex=2)
par(cex.axis=1.2)


  plot(c(200, 1500), c(0, 1), t = "n", 
       ylab = "cumsum of % deaths in window (relative to population)", xlab = "Lifespan Cut-off Age (days)", yaxt="n", xaxs="i", main = "Effect Curve")
  msequence <- seq(20, 1500, 15)
  std <- function(x) sd(x)/sqrt(length(x))
  mm <- c()
  for(day in msequence){
    day2 <- day + 15
    nCdR <- length(which(C > day & C <= day2)) / length(c(C,D))
    nDdR <- length(which(D > day & D <= day2)) / length(c(C,D))
    mm <- rbind(mm, c(nCdR, nDdR))
  }
  abline(h = seq(0, 15, 1), lty=2, col = "gray")
  points(msequence, cumsum(mm[,1]), t = 'l', col = "#FF3399", lwd=2)
  points(msequence, cumsum(mm[,2]), t = 'l', col = "#CCCC99", lwd=2)
  #polygon(c(msequence, rev(msequence)), c(mm[,1] + mm[,2], rev(mm[,1] - mm[,2])), col = rgb(0,0,1,0.5), border = NA)

  #points(msequence, mm[,3], t = 'l', col = rgb(.4,0.2,1,1), lwd=2)
  #polygon(c(msequence, rev(msequence)), c(mm[,3] + mm[,4], rev(mm[,3] - mm[,4])), col = rgb(.4,0.2,1,0.5), border = NA)


  legend("topleft", c("H", "D"), col = c("#FF3399", "#CCCC99"), pch=18)
  axis(2, at = seq(0, 1, 0.1), seq(0, 1, 0.1), las=2)
  #dev.off()

  abline(v = 970)
  abline(v = 520)
dev.off()


hist(mm[,1] - mm[,2], col = rgb(1,0,0,0.5), breaks = seq(0, 0.02, .0001))
hist(mm[,2], add = TRUE, col = rgb(0,1,0,0.5), breaks = seq(0, 0.02, .0001))


setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny")
pdf("Vita2B_RelativeDeaths_15day_Population.pdf", width = 16, height = 12)

par(cex=2)
par(cex.axis=1.2)


plot(msequence, 1000*smooth(mm[,1]-mm[,2]), t = "b", pch=19, 
     col= c("#FF3399", "#CCCC99")[1-as.numeric(smooth(mm[,1]-mm[,2]) > 0)+1], xaxt="n", ylab = "Difference in % deaths C versus H in windows", main = "Vita2b")
axis(1, at = seq(0, 1600, 200), seq(0, 1600, 200), las=1)
abline(h = 0)
legend("topleft", c("H", "D"), col = c("pink", "brown"), pch=18)

dev.off()

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny")
pdf("Vita2B_Deaths_15day_OwnAllele.pdf", width = 16, height = 12)

par(cex=2)
par(cex.axis=1.2)


  plot(c(200, 1500), c(0, 1), t = "n", ylab = "cumsum of % deaths in window (relative to same allele)", xlab = "Lifespan Cut-off Age (days)", yaxt="n", xaxs="i", main = "Effect Curve")
  msequence <- seq(20, 1500, 15)
  std <- function(x) sd(x)/sqrt(length(x))
  mm <- c()
  for(day in msequence){
    day2 <- day + 15
    nCdR <- length(which(C > day & C <= day2)) / length(c(C))
    nDdR <- length(which(D > day & D <= day2)) / length(c(D))
    mm <- rbind(mm, c(nCdR, nDdR))
  }
  abline(h = seq(0, 15, 1), lty=2, col = "gray")
  points(msequence, cumsum(mm[,1]), t = 'l', col = "#FF3399", lwd=2)
  points(msequence, cumsum(mm[,2]), t = 'l', col = "#CCCC99", lwd=2)
  #polygon(c(msequence, rev(msequence)), c(mm[,1] + mm[,2], rev(mm[,1] - mm[,2])), col = rgb(0,0,1,0.5), border = NA)

  #points(msequence, mm[,3], t = 'l', col = rgb(.4,0.2,1,1), lwd=2)
  #polygon(c(msequence, rev(msequence)), c(mm[,3] + mm[,4], rev(mm[,3] - mm[,4])), col = rgb(.4,0.2,1,0.5), border = NA)


  legend("topleft", c("H", "D"), col = c("#FF3399", "#CCCC99"), pch=18)
  axis(2, at = seq(0, 1, 0.1), seq(0, 1, 0.1), las=2)
  #dev.off()

  abline(v = 970)
  abline(v = 520)
dev.off()


