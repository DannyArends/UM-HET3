setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

bw <- read.table("42day_bodyweight.txt",sep="\t",na.strings=c("", "NA", "x"), header=TRUE, row.names=1)

gtsp <- pull.genoprob(mcross)
cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    adjLongevity = NA, 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw6 = bw[pull.pheno(mcross)[, "GenoID"],"BW_42d"],
                    adjBw6 = NA,
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
idx <- which(cdata[, "longevity"] >= (365/2) & !is.na(cdata[, "bw6"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw6 <- lm(bw6 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw6"] <- round(as.numeric(coef(lm.null.bw6)["(Intercept)"]) + residuals(lm.null.bw6), 2)


all <- c("4_154254581", "11_5628810")
names(all) <- c("Bwle4b", "Bwle11a")

genotypes <- c()
for(marker in all){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  genotypes <- cbind(genotypes, gts)
}
colnames(genotypes) <- all

### TODO add the error bars of Rho

doPoints <- function(gts, stepsize = 2, col = 1, sex = c(0,1)) {
  for(ii in seq(12, 60, stepsize)){
    iR <- which(cdata[gts, "adjBw6"] > ii & cdata[gts, "adjBw6"] < ii + stepsize & cdata[gts,"sex"] %in% sex)
    mm <- mean(cdata[gts[iR], "adjLongevity"])
    ms <- sd(cdata[gts[iR], "adjLongevity"]) / sqrt(length(iR))
    p <- ii + (stepsize/2) + seq(-0.3,0.3,0.15)[-3][col]
    if(length(iR) > 2){
      s <- round(log2(length(iR))) - 1
      points(p, mm, pch = 19, col = col.main[col], cex=s)
      points(c(p,p), c(mm - ms, mm + ms), pch = 19, col = col.main[col], t = "l", lty=2)
    }
  }

  iR <- which(cdata[gts,"sex"] %in% sex)
  abline(lm(cdata[gts[iR], "adjLongevity"] ~ cdata[gts[iR], "adjBw6"]), col = col.main[col], lwd=4, untf=T)
}

for(i in 1:length(all)){
  x <- all[i]
  CH <- which(genotypes[,x] == "AC")
  BH <- which(genotypes[,x] == "BC")
  CD <- which(genotypes[,x] == "AD")
  BD <- which(genotypes[,x] == "BD")

  col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
  col.main <- adjustcolor( col.main, alpha.f = 0.6)

  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")
  pdf(paste0("BW42_CTL_",names(all)[i],".pdf"), width = 36, height = 12)
  op <- par(mfrow = c(1,3))
  par(cex=1.5)
  par(cex.axis=1.2)

    plot(c(12, 26), c(500, 1100), t = "n", xlab = "Adjusted bodyweight (42d)", ylab = "Adjusted longevity", 
         main = paste0("CTL [All] @ ", x), log = "y", yaxt = "n")

    doPoints(CH, 2, 1)
    doPoints(CD, 2, 2)
    doPoints(BH, 2, 3)
    doPoints(BD, 2, 4)

    axis(2, at = seq(500, 1100, 100), seq(500, 1100, 100), las=2)
    legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=19, bg = "white", ncol=4, bty = "n")

    plot(c(12, 26), c(500, 1100), t = "n", xlab = "Adjusted bodyweight (42d)", ylab = "Adjusted longevity", 
         main = paste0("CTL [Males] @ ", x), log = "y", yaxt = "n")

    doPoints(CH, 2, 1,1)
    doPoints(CD, 2, 2,1)
    doPoints(BH, 2, 3,1)
    doPoints(BD, 2, 4,1)

    axis(2, at = seq(500, 1100, 100), seq(500, 1100, 100), las=2)
    legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=19, bg = "white", ncol=4, bty = "n")

    plot(c(12, 26), c(500, 1100), t = "n", xlab = "Adjusted bodyweight (42d)", ylab = "Adjusted longevity", 
         main = paste0("CTL [Females] @ ", x), log = "y", yaxt = "n")

    doPoints(CH, 2, 1, 0)
    doPoints(CD, 2, 2, 0)
    doPoints(BH, 2, 3, 0)
    doPoints(BD, 2, 4, 0)

    axis(2, at = seq(500, 1100, 100), seq(500, 1100, 100), las=2)
    legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=19, bg = "white", ncol=4, bty = "n")

  dev.off()
}

doPoints <- function(gts, stepsize = 2, col = 1, sex = c(0,1)) {
  for(ii in seq(20, 60, stepsize)){
    iR <- which(cdata[gts, "adjBw6"] > ii & cdata[gts, "adjBw6"] < ii + stepsize & cdata[gts,"sex"] %in% sex)
    mm <- mean(cdata[gts[iR], "adjLongevity"])
    ms <- sd(cdata[gts[iR], "adjLongevity"]) / sqrt(length(iR))
    p <- ii + (stepsize/2) + seq(-0.3,0.3,0.15)[-3][col]
    if(length(iR) > 2){
      s <- round(log2(length(iR))) - 1
      points(p, mm, pch = 19, col = col.main[col], cex=s)
      points(c(p,p), c(mm - ms, mm + ms), pch = 19, col = col.main[col], t = "l", lty=2)
    }
  }

  iR <- which(cdata[gts,"sex"] %in% sex)
  abline(lm(cdata[gts[iR], "adjLongevity"] ~ cdata[gts[iR], "adjBw6"]), col = col.main[col], lwd=4, untf=T)
}


for(i in 1:length(all)){
  x <- all[i]
  CH <- which(genotypes[,x] == "AC")
  BH <- which(genotypes[,x] == "BC")
  CD <- which(genotypes[,x] == "AD")
  BD <- which(genotypes[,x] == "BD")

  col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
  col.main <- adjustcolor( col.main, alpha.f = 0.6)

  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")
  pdf(paste0("Figure_3_CTL_",names(all)[i],"_bin3_v2.pdf"), width = 36, height = 12)
  op <- par(mfrow = c(1,3))
  par(cex=1.5)
  par(cex.axis=1.2)

    stepsize = 3

    plot(c(20, 48), c(500, 1100), t = "n", xlab = "Adjusted bodyweight (6mo)", ylab = "Adjusted longevity", 
         main = paste0("CTL [All] @ ", x), yaxt = "n")

    doPoints(CH, stepsize, 1)
    doPoints(CD, stepsize, 2)
    doPoints(BH, stepsize, 3)
    doPoints(BD, stepsize, 4)

    axis(2, at = seq(500, 1100, 100), seq(500, 1100, 100), las=2)
    legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=19, bg = "white", ncol=4, bty = "n")

    plot(c(20, 48), c(500, 1100), t = "n", xlab = "Adjusted bodyweight (6mo)", ylab = "Adjusted longevity", 
         main = paste0("CTL [Males] @ ", x), yaxt = "n")

    doPoints(CH, stepsize, 1,1)
    doPoints(CD, stepsize, 2,1)
    doPoints(BH, stepsize, 3,1)
    doPoints(BD, stepsize, 4,1)

    axis(2, at = seq(500, 1100, 100), seq(500, 1100, 100), las=2)
    legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=19, bg = "white", ncol=4, bty = "n")

    plot(c(20, 48), c(500, 1100), t = "n", xlab = "Adjusted bodyweight (6mo)", ylab = "Adjusted longevity", 
         main = paste0("CTL [Females] @ ", x), yaxt = "n")

    doPoints(CH, stepsize, 1, 0)
    doPoints(CD, stepsize, 2, 0)
    doPoints(BH, stepsize, 3, 0)
    doPoints(BD, stepsize, 4, 0)

    axis(2, at = seq(500, 1100, 100), seq(500, 1100, 100), las=2)
    legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=19, bg = "white", ncol=4, bty = "n")

  dev.off()
}




