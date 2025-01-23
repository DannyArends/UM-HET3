setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

gtsp <- pull.genoprob(mcross)
cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    adjLongevity = NA, 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw6 = as.numeric(pull.pheno(mcross)[, "bw6"]),
                    adjBw6 = NA,
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
idx <- which(cdata[, "longevity"] >= 365 & !is.na(cdata[, "bw6"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw6 <- lm(bw6 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw6"] <- round(as.numeric(coef(lm.null.bw6)["(Intercept)"]) + residuals(lm.null.bw6), 2)


all <- c("1_3010274", "2_65326540", "3_159581164", "4_11234654", "5_34118900", "6_128506813",  "7_16072018", "8_95039510", "10_18144599", "13_53167285", "13_76135291", "15_3288506", "15_79892499", "15_99306167", "17_35023240")
names(all) <- c("Soma1a", "Soma2a", "Soma3a", "Soma4a", "Soma5a", "Soma6a", "Soma7a", "Soma8a", "Soma10a", "Soma13a", "Soma13b", "Soma15a", "Soma15b", "Soma15c", "Soma17a")

genotypes <- c()
for(marker in all){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  genotypes <- cbind(genotypes, gts)
}
colnames(genotypes) <- all

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


for(i in length(all)[1]){
  x <- all[i]
  CH <- which(genotypes[,x] == "AC")
  BH <- which(genotypes[,x] == "BC")
  CD <- which(genotypes[,x] == "AD")
  BD <- which(genotypes[,x] == "BD")


  col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27") # green, blue, purple, red
  col.main <- adjustcolor( col.main, alpha.f = 0.6)

  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")
  pdf(paste0("SOMA_",names(all)[i],"_SCALE.pdf"), width = 36, height = 12)
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
    legend("topright", c("CH", "CD", "BH", "BD")[c(1, 2, 3, 4)], col = col.main, pch=19, bg = "white", ncol=4, bty = "n")

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
    legend("topleft", c("3-5", "6-11", "12-22", "23-45", "46-90", "91-181", "182-362"), col = "black", pch=19, bg = "white", ncol=1, bty = "n", pt.cex=1:7)
  dev.off()
}




