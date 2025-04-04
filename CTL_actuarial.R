#
# Make all of them for 42, 12mo, 18mo & 24mo
#
#
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
idx <- which(cdata[, "longevity"] >= 0 & !is.na(cdata[, "bw6"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw6 <- lm(bw6 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw6"] <- round(as.numeric(coef(lm.null.bw6)["(Intercept)"]) + residuals(lm.null.bw6), 2)

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

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny/SOMA Effects_185_day")

all <- c("1_3010272","1_86216552","2_13600088", "2_60201233","2_161871392","3_87974845","3_159581164","4_30761996","4_107374161","6_8006720",
         "6_138658041","7_16072018","7_120086292","8_71684276","8_126505019","9_51116640","10_18144599","11_97448477","12_71677220",
         "12_118179607","13_19367506","13_98521647","14_30957748", "14_101437466","15_3288506","16_75758401","17_26542857",
         "18_52488251","19_3403302","19_53851357")
names(all) <- c("Soma1a","Soma1b","Soma2a","Soma2b","Soma2c","Soma3a","Soma3b","Soma4a","Soma4b","Soma6a","Soma6b","Soma7a","Soma7b","Soma8a",
                "Soma8b","Soma9a","Soma10a","Soma11a","Soma12a","Soma12b","Soma13a","Soma13b","Soma14a","Soma14b","Soma15a","Soma16a",
                "Soma17a","Soma18a","Soma19a","Soma19b")

for(m in names(all)){
  pdf(paste0(m, "_BW185_actuarial.pdf"), width = (28+14) * .8, height = 12)
  op <- par(mfrow=c(1,3))
  op <- par(cex = 2)
  par(mar = c(5, 4, 4, 2))
  for(s in list(c(0,1), 0, 1)){
    selected <- all[m]
    gts <- genotypes[, selected]
    sex <- s

    corM <- c()
    allN <- c()
    confL <- c()
    confU <- c()
    for(x in seq(185, 1100, 15)){
      CH <- which(gts == "AC" & cdata[, "sex"] %in% sex & cdata[, "longevity"] >= x)
      BH <- which(gts == "BC" & cdata[, "sex"] %in% sex & cdata[, "longevity"] >= x)
      CD <- which(gts == "AD" & cdata[, "sex"] %in% sex & cdata[, "longevity"] >= x)
      BD <- which(gts == "BD" & cdata[, "sex"] %in% sex & cdata[, "longevity"] >= x)
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

  #  col.main <- c("#01A654", "#1750A3", "#714F99", "#EE3129")
    col.main <- c("#00A654", "#B16BE6", "#004BAD", "#F02D27") # green, blue, purple, red
    add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
    col.alpha2 <- add.alpha(col.main, 0.1)

    ylim <- -0.45
    if(length(sex) == 1 && sex == 0) ylim <- -0.45

    pop <- " Combined"
    if(length(sex) == 1 && sex == 0) pop <- " Female"
    if(length(sex) == 1 && sex == 1) pop <- " Male"

#### labels = bquote(atop(bold(bolditalic(.(m1)) ~ x ~ bolditalic(.(m2))), bold(at ~ T[.(timepoint)])))

    plot(c(185, 1000), c(ylim, 0.2), t = "n", xlab = "", ylab = "", 
         main = "", yaxt = "n", yaxs = "i")
    title(bquote(atop(bold(bolditalic(.(m)) ~ .(pop)))), line = 0.2, cex = 1.5)
    title(ylab=bquote(atop(bold("BW185 vs. T-age Correlation"))), line=1.75)
    title(xlab=bquote(atop(bold("Truncation age [d]"))), line=3.5)

    #abline(v = 365, col = "gray")
    #abline(h = -0.3, col = "gray")
    #abline(h = -0.25, col = "gray")
    #abline(h = -0.15, col = "gray")
    text(x=900, y = ylim+0.1, labels = paste0("max LOD = ", max(lods)), cex=.8)

    points(seq(185, 1100, 15), corM[,1], t = "l", col = col.main[1], lwd = 2)
    polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,1], rev(confU[,1])),col = col.alpha2[1], border=NA)

    points(seq(185, 1100, 15), corM[,2], t = "l", col = col.main[2], lwd = 2)
    polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,2], rev(confU[,2])),col = col.alpha2[2], border=NA)

    points(seq(185, 1100, 15), corM[,3], t = "l", col = col.main[3], lwd = 2)
    polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,3], rev(confU[,3])),col = col.alpha2[3], border=NA)

    points(seq(185, 1100, 15), corM[,4], t = "l", col = col.main[4], lwd = 2)
    polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,4], rev(confU[,4])),col = col.alpha2[4], border=NA)
    axis(2, at = seq(-.5, .2, .1), las=2)
    axis(2, at = seq(-.5, .2, .05), rep("", length(seq(-.5, .2, .05))), las=2)
    axis(1, at = seq(100, 1100, 100), rep("",length(seq(100, 1100, 100))))
    abline(h = ylim + (2.75/50), lty=1, col = "lightgray")
    points(seq(185, 1100, 15), (lods / 50) + ylim, t = "l", lwd=2)
    legend("topleft", c("CH", "BH", "CD", "BD")[c(1,3,2,4)], col = col.main[c(1,3,2,4)], lwd=4, bg = "white", ncol=4, bty = "n")
    axis(4, at = c(ylim + 0.04, ylim + 0.08), c(2,4), las = 2,  cex.axis =1.0)
  }
  dev.off()
}

##significance











