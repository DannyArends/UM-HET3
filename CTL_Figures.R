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
idx <- which(cdata[, "longevity"] >= 185 & !is.na(cdata[, "bw6"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw6 <- lm(bw6 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw6"] <- round(as.numeric(coef(lm.null.bw6)["(Intercept)"]) + residuals(lm.null.bw6), 2)

## TODO 

all <- c("1_3010272","1_86216552","2_13600088", "2_60201233","2_161871392","3_87974845","3_159581164","4_30761996","4_107374161","6_8006720",
         "6_138658041","7_16072018","7_120086292","8_71684276","8_126505019","9_51116640","10_18144599","11_97448477","12_71677220",
         "12_118179607","13_19367506","13_98521647","14_30957748", "14_101437466","15_3288506","16_75758401","17_26542857",
         "18_52488251","19_3403302","19_53851357")
names(all) <- c("Soma1a","Soma1b","Soma2a","Soma2b","Soma2c","Soma3a","Soma3b","Soma4a","Soma4b","Soma6a","Soma6b","Soma7a","Soma7b","Soma8a",
                "Soma8b","Soma9a","Soma10a","Soma11a","Soma12a","Soma12b","Soma13a","Soma13b","Soma14a","Soma14b","Soma15a","Soma16a",
                "Soma17a","Soma18a","Soma19a","Soma19b")

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

  #setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")
  #pdf(paste0("SOMA_",names(all)[i],"_SCALE.pdf"), width = 36, height = 12)
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
  #dev.off()
}




