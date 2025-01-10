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
idx <- which(cdata[, "longevity"] >= (365/2) & !is.na(cdata[, "bw6"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw6 <- lm(bw6 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw6"] <- round(as.numeric(coef(lm.null.bw6)["(Intercept)"]) + residuals(lm.null.bw6), 2)

all <- c("1_3010274", "2_65326540", "3_159581164", "4_11234654", "5_34118900", "6_128506813",  "7_16072018", "8_95039510", "10_18144599", "13_53167285", "13_76135291", "15_3288506", "15_79892499", "15_99306167", "17_35023240")
names(all) <- c("Bwle1a", "Bwle2a", "Bwle3a", "Bwle4a", "Bwle5a", "Bwle6a", "Bwle7a", "Bwle8a", "Bwle10a", "Bwle13a", "Bwle13b", "Bwle15a", "Bwle15b", "Bwle15c", "Bwle17a")

genotypes <- c()
for(marker in all){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  genotypes <- cbind(genotypes, gts)
}
colnames(genotypes) <- all

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")


toP <- function(allCor, allN){
  ## correlation differences to P-value / LOD scores
  pC <- c()
  for(x in 1:nrow(allCor)){
    pvs <- c()

    # Compute the allele deltas relative to eachother
    cor <- na.omit(allCor[x,])
    z <- .5*log((1.0 + cor)/(1.0 - cor))
    n <- allN[x,]
    df <- n-(length(cor)-1)
    sumOfSq <- sum(df * z^2)
    sqOfSum <- sum(df * z)
    denom <- sum(df)
    Cv <- sumOfSq - (sqOfSum^2/ denom)
    pC <- c(pC, pchisq(Cv, 1, 0, FALSE))
  }
  names(pC) <- rownames(allCor)
  return(list(pC))
}



for(m in names(all)){
  pdf(paste0(m, "_BW185_actuarial.pdf"), width = 28+14, height = 12)
  op <- par(mfrow=c(1,3))
  op <- par(cex = 2)
  par(mar = c(5, 5, 4, 3))
  for(s in list(c(0,1), 1, 0)){
    selected <- all[m]
    gts <- genotypes[, selected]
    sex <- s

    corM <- c()
    allN <- c()
    confL <- c()
    confU <- c()
    for(x in seq(185, 1100, 15)){
      CH <- which(gts == "AC" & cdata[, "sex"] %in% sex & cdata[, "adjLongevity"] > x)
      BH <- which(gts == "BC" & cdata[, "sex"] %in% sex & cdata[, "adjLongevity"] > x)
      CD <- which(gts == "AD" & cdata[, "sex"] %in% sex & cdata[, "adjLongevity"] > x)
      BD <- which(gts == "BD" & cdata[, "sex"] %in% sex & cdata[, "adjLongevity"] > x)
      cCH <- NA; cBH <- NA; cCD <- NA; cBD <- NA

      if(length(CH) > 100){
        cCH <- cor(cdata[CH, "adjLongevity"], cdata[CH, "adjBw6"], use = "pair", method = "spearman");
        confCH <- cor.test(cdata[CH, "adjLongevity"], cdata[CH, "adjBw6"], use = "pair", method = "pearson", conf.level = 0.50);
      }
      if(length(BH) > 100){
        cBH <- cor(cdata[BH, "adjLongevity"], cdata[BH, "adjBw6"], use = "pair", method = "spearman");
        confBH <- cor.test(cdata[BH, "adjLongevity"], cdata[BH, "adjBw6"], use = "pair", method = "pearson", conf.level = 0.50);
      }
      if(length(CD) > 100){
        cCD <- cor(cdata[CD, "adjLongevity"], cdata[CD, "adjBw6"], use = "pair", method = "spearman");
        confCD <- cor.test(cdata[CD, "adjLongevity"], cdata[CD, "adjBw6"], use = "pair", method = "pearson", conf.level = 0.50);
      }
      if(length(BD) > 100){
        cBD <- cor(cdata[BD, "adjLongevity"], cdata[BD, "adjBw6"], use = "pair", method = "spearman");
        confBD <- cor.test(cdata[BD, "adjLongevity"], cdata[BD, "adjBw6"], use = "pair", method = "pearson", conf.level = 0.50);
      }

      allN <- rbind(allN, c(length(CH), length(BH), length(CD), length(BD)))
      corM <- rbind(corM, c(cCH,cBH, cCD, cBD))

      confL <- rbind(confL, c(round(as.numeric(unlist(confCH)["conf.int1"]),2),
                              round(as.numeric(unlist(confBH)["conf.int1"]),2), 
                              round(as.numeric(unlist(confCD)["conf.int1"]),2), 
                              round(as.numeric(unlist(confBD)["conf.int1"]),2)))
      confU <- rbind(confU, c(round(as.numeric(unlist(confCH)["conf.int2"]),2),
                              round(as.numeric(unlist(confBH)["conf.int2"]),2), 
                              round(as.numeric(unlist(confCD)["conf.int2"]),2), 
                              round(as.numeric(unlist(confBD)["conf.int2"]),2)))

      for(x in 1:nrow(corM)){
        for(y in 1:ncol(corM)){
          mid <- corM[x, y]
          adj <- (confU[x,y] - confL[x, y]) / 2

          confU[x,y] <- mid + adj
          if(is.na(confU[x,y])) confU[x,y] <- 0
          confL[x,y] <- mid - adj
          if(is.na(confL[x,y])) confL[x,y] <- 0
        }
      }
    }

    ttt <-toP(corM, allN)
    lods <- round(-log10(ttt[[1]]),2)

    col.main <- c("#01A654", "#1750A3", "#714F99", "#EE3129")
    add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
    col.alpha <- add.alpha(col.main, 0.1)
    col.alpha2 <- add.alpha(col.main, 0.1)

    ylim <- -0.5
    if(length(sex) == 1 && sex == 0) ylim <- -0.35

    pop <- "Combined"
    if(length(sex) == 1 && sex == 0) pop <- "Female"
    if(length(sex) == 1 && sex == 1) pop <- "Male"
    main <- paste0(pop, " at ", m)

    plot(c(185, 1100), c(ylim, 0.1), t = "n", xlab = "Truncation age (days)", ylab = "Correlation BW185 to Tage", 
         main = main,yaxt = "n", yaxs = "i")

    text(x=1000, y = ylim+0.1, labels = paste0("max LOD = ", max(lods)), cex=.8)

    points(seq(185, 1100, 15), corM[,1], t = "l", col = col.main[1], lwd = 2)
    polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,1], rev(confU[,1])),col = col.alpha2[1], border=NA)

    points(seq(185, 1100, 15), corM[,2], t = "l", col = col.main[2], lwd = 2)
    polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,2], rev(confU[,2])),col = col.alpha2[2], border=NA)

    points(seq(185, 1100, 15), corM[,3], t = "l", col = col.main[3], lwd = 2)
    polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,3], rev(confU[,3])),col = col.alpha2[3], border=NA)

    points(seq(185, 1100, 15), corM[,4], t = "l", col = col.main[4], lwd = 2)
    polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,4], rev(confU[,4])),col = col.alpha2[4], border=NA)
    axis(2, at = seq(-.5, .1, .1), las=2)
    axis(2, at = seq(-.5, .1, .05), rep("", length(seq(-.5, .1, .05))), las=2)
    axis(1, at = seq(100, 1100, 100), rep("",length(seq(100, 1100, 100))))
    abline(h = ylim + (2.75/50), lty=1, col = "lightgray")
    points(seq(185, 1100, 15), (lods / 50) + ylim, t = "l", lwd=2)
    legend("topleft", c("CH", "CD", "BH", "BD"), col = col.main, lwd=4, bg = "white", ncol=4, bty = "n")
    axis(4, at = c(ylim + 0.04, ylim + 0.08), c(2,4), las = 2,  cex.axis =1.0)
  }
  dev.off()
}

##significance











