#
# CTL_Figures_v1.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Create actuarial CTL plots per Soma locus across the different time points (Supplemental files)
#

library(qtl)

source("ActuarialMapping/adjustXprobs.R")
mcross <- read.cross(format="csvr", file="DataSet/um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

adjustPHE <- function(cdata, days = 0, column = "bw42", out = c("is42", "adjBw42", "adjLs42")){
  iz <- cdata[, "longevity"] >= days & !is.na(cdata[, column])
  mdata <- cdata[which(iz),]

  cdata[, out[1]] <<- iz

  lm.null <- lm(as.formula(paste0(column, " ~ sex + site + cohort + treatment")), data = mdata)
  cdata[which(iz), out[2]] <<- round(as.numeric(coef(lm.null)["(Intercept)"]) + residuals(lm.null), 2)

  lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = mdata)
  cdata[which(iz), out[3]] <<- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)
}

#42 days
bwA <- read.csv("DataSet/bodyweights/bw_RichM.txt", sep = "\t", header = TRUE, comment.char = "#", row.names=1, na.strings = c("NA", "", "x"))
bw <- read.csv("DataSet/bodyweights/ITP_50601.csv", header = TRUE, comment.char = "#", skip=11, row.names=2,na.strings = c("NA", "", "x"))
bw <- bw[which(bw[, "DA2024"] == 1),]

trait <- read.csv("DataSet/bodyweights/ITP_10003.csv", header = TRUE, comment.char = "#", skip=10, row.names=2,na.strings = c("NA", "", "x"))
trait <- trait[which(trait[, "DA2024"] == 1),]
snames <- as.character(pull.pheno(mcross)[, "GenoID"])

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw42 = as.numeric(bw[snames, "Value"]),
                    bw6 = as.numeric(pull.pheno(mcross)[, "bw6"]),
                    bw12 = as.numeric(trait[snames, "Value"]),
                    bw18 = as.numeric(pull.pheno(mcross)[, "bw18"]),
                    bw24 = as.numeric(pull.pheno(mcross)[, "bw24"]),
                    is42 = NA,
                    adjBw42 = NA,
                    adjLs42 = NA, 

                    is6 = NA,
                    adjBw6 = NA,
                    adjLs6 = NA, 

                    is12 = NA,
                    adjBw12 = NA,
                    adjLs12 = NA, 

                    is18 = NA,
                    adjBw18 = NA,
                    adjLs18 = NA, 

                    is24 = NA,
                    adjBw24 = NA,
                    adjLs24 = NA, 
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

rownames(cdata) <- snames

noData <- snames[which(is.na(cdata[, "bw42"]))]
addData <- noData[which(noData %in% rownames(bwA))]

#Add the data from Rich to the data we had
cdata[addData, "bw42"] <- bwA[addData, "W42d"]

adjustPHE(cdata, 0, "bw42", c("is42", "adjBw42", "adjLs42"))
adjustPHE(cdata, 185, "bw6", c("is6", "adjBw6", "adjLs6"))
adjustPHE(cdata, 365, "bw12", c("is12", "adjBw12", "adjLs12"))
adjustPHE(cdata, 550, "bw18", c("is18", "adjBw18", "adjLs18"))
adjustPHE(cdata, 725, "bw24", c("is24", "adjBw24", "adjLs24"))

### All our 30 Vita Loci
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

timepoints <- c(42, 185, 365, 550, 730)
names(timepoints) <- c("42", "6", "12", "18", "24")

for(tp in c("42", "6", "12", "18", "24")){
  longCol <- paste0("adjLs", tp)
  bwCol <- paste0("adjBw", tp)
  bwStCol <- paste0("bw", tp)
  startTP <- timepoints[tp]
  for(m in names(all)){
    pdf(paste0("DataSet/output/", m, "_", tp,"_actuarial.pdf"), width = (28+14) * .8, height = 12)
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
      for(x in seq(startTP, 1100, 15)){
        CH <- which(gts == "AC" & cdata[, "sex"] %in% sex & cdata[, "longevity"] >= x & !is.na(cdata[, bwStCol]))
        BH <- which(gts == "BC" & cdata[, "sex"] %in% sex & cdata[, "longevity"] >= x & !is.na(cdata[, bwStCol]))
        CD <- which(gts == "AD" & cdata[, "sex"] %in% sex & cdata[, "longevity"] >= x & !is.na(cdata[, bwStCol]))
        BD <- which(gts == "BD" & cdata[, "sex"] %in% sex & cdata[, "longevity"] >= x & !is.na(cdata[, bwStCol]))
        cCH <- NA; cBH <- NA; cCD <- NA; cBD <- NA

        if(length(CH) > 100){
          cCH <- cor(cdata[CH, longCol], cdata[CH, bwCol], use = "pair", method = "spearman");
          confCH <- cor.test(cdata[CH, longCol], cdata[CH, bwCol], use = "pair", method = "pearson", conf.level = 0.50);
        }
        if(length(BH) > 100){
          cBH <- cor(cdata[BH, longCol], cdata[BH, bwCol], use = "pair", method = "spearman");
          confBH <- cor.test(cdata[BH, longCol], cdata[BH, bwCol], use = "pair", method = "pearson", conf.level = 0.50);
        }
        if(length(CD) > 100){
          cCD <- cor(cdata[CD, longCol], cdata[CD, bwCol], use = "pair", method = "spearman");
          confCD <- cor.test(cdata[CD, longCol], cdata[CD, bwCol], use = "pair", method = "pearson", conf.level = 0.50);
        }
        if(length(BD) > 100){
          cBD <- cor(cdata[BD, longCol], cdata[BD, bwCol], use = "pair", method = "spearman");
          confBD <- cor.test(cdata[BD, longCol], cdata[BD, bwCol], use = "pair", method = "pearson", conf.level = 0.50);
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

      plot(c(42, 1000), c(ylim, 0.2), t = "n", xlab = "", ylab = "", 
           main = "", yaxt = "n", yaxs = "i")
      title(bquote(atop(bold(bolditalic(.(m)) ~ .(pop)) ~ "(" ~ .(startTP) ~ d ~ Mass ~ ")")), line = 0.2, cex = 1.5)
      title(ylab=bquote(atop(bold("BW185 vs. T-age Correlation"))), line=1.75)
      title(xlab=bquote(atop(bold("Truncation age [d]"))), line=3.5)

      #abline(v = 365, col = "gray")
      #abline(h = -0.3, col = "gray")
      #abline(h = -0.25, col = "gray")
      #abline(h = -0.15, col = "gray")
      text(x=900, y = ylim+0.1, labels = paste0("max -logP = ", max(lods)), cex=.8)

      points(seq(startTP, 1100, 15), corM[,1], t = "l", col = col.main[1], lwd = 2)
      polygon(c(seq(startTP, 1100, 15), rev(seq(startTP, 1100, 15))), c(confL[,1], rev(confU[,1])),col = col.alpha2[1], border=NA)

      points(seq(startTP, 1100, 15), corM[,2], t = "l", col = col.main[2], lwd = 2)
      polygon(c(seq(startTP, 1100, 15), rev(seq(startTP, 1100, 15))), c(confL[,2], rev(confU[,2])),col = col.alpha2[2], border=NA)

      points(seq(startTP, 1100, 15), corM[,3], t = "l", col = col.main[3], lwd = 2)
      polygon(c(seq(startTP, 1100, 15), rev(seq(startTP, 1100, 15))), c(confL[,3], rev(confU[,3])),col = col.alpha2[3], border=NA)

      points(seq(startTP, 1100, 15), corM[,4], t = "l", col = col.main[4], lwd = 2)
      polygon(c(seq(startTP, 1100, 15), rev(seq(startTP, 1100, 15))), c(confL[,4], rev(confU[,4])),col = col.alpha2[4], border=NA)
      axis(2, at = seq(-.5, .2, .1), las=2)
      axis(2, at = seq(-.5, .2, .05), rep("", length(seq(-.5, .2, .05))), las=2)
      axis(1, at = seq(100, 1100, 100), rep("",length(seq(100, 1100, 100))))
      abline(h = ylim + (2.75/50), lty=1, col = "lightgray")
      points(seq(startTP, 1100, 15), (lods / 50) + ylim, t = "l", lwd=2)
      legend("topleft", c("CH", "BH", "CD", "BD")[c(1,3,2,4)], col = col.main[c(1,3,2,4)], lwd=4, bg = "white", ncol=4, bty = "n")
      axis(4, at = c(ylim + 0.04, ylim + 0.08), c(2,4), las = 2,  cex.axis =1.0)
    }
    dev.off()
  }
}

