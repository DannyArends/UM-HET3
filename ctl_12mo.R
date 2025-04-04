setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

#
# Use a different clor scheme (blue - red)
#

### TODO: CTL for Males and Females separately (Done)

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

trait <- read.csv("ITP_10003.csv", header = TRUE, comment.char = "#", skip=10, row.names=2,na.strings = c("NA", "", "x"))
trait <- trait[which(trait[, "DA2024"] == 1),]
snames <- as.character(pull.pheno(mcross)[, "GenoID"])


gtsp <- pull.genoprob(mcross)
cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    adjLongevity = NA, 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw12 = as.numeric(trait[snames, "Value"]),
                    adjBw12 = NA,
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
rownames(cdata) <- pull.pheno(mcross)[, "GenoID"]

idx <- which(cdata[, "longevity"] >= 365 & !is.na(cdata[, "bw12"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw12 <- lm(bw12 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw12"] <- round(as.numeric(coef(lm.null.bw12)["(Intercept)"]) + residuals(lm.null.bw12), 2)


### Plot males and females
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

corM <- c()
allN <- c()
confL <- c()
confU <- c()
for(x in seq(365, 1100, 15)){
  male <- which(cdata[, "sex"] %in% 0 & cdata[, "adjLongevity"] > x)
  fema <- which(cdata[, "sex"] %in% 1 & cdata[, "adjLongevity"] > x)
  cMale <- NA; cFema <- NA;

  if(length(male) > 100){
    cMale <- cor(cdata[male, "adjLongevity"], cdata[male, "adjBw12"], use = "pair", method = "spearman");
    confMale <- cor.test(cdata[male, "adjLongevity"], cdata[male, "adjBw12"], use = "pair", method = "pearson", conf.level = 0.50);
  }
  if(length(fema) > 100){
    cFema <- cor(cdata[fema, "adjLongevity"], cdata[fema, "adjBw12"], use = "pair", method = "spearman");
    confFema <- cor.test(cdata[fema, "adjLongevity"], cdata[fema, "adjBw12"], use = "pair", method = "pearson", conf.level = 0.50);
  }
  allN <- rbind(allN, c(length(male), length(fema)))
  corM <- rbind(corM, c(cMale, cFema))

  confL <- rbind(confL, c(round(as.numeric(unlist(confMale)["conf.int1"]),2),
                          round(as.numeric(unlist(confFema)["conf.int1"]),2)))
  confU <- rbind(confU, c(round(as.numeric(unlist(confMale)["conf.int2"]),2),
                          round(as.numeric(unlist(confFema)["conf.int2"]),2)))
  # Adjust Conf
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

col.main <- c("#FF3333", "#00AEEF")
add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
col.alpha2 <- add.alpha(col.main, 0.1)

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")
pdf(paste0("CTL_MF_12months.pdf"), width = 14, height = 12)
op <- par(cex = 2)

  plot(c(42, 1100), c(-0.5, 0.2), t = "n", xlab = "Truncation age (days)", ylab = "Correlation BW365 to Tage", 
       main = "Correlation Male/Female Bodyweight 12-months",yaxt = "n", yaxs = "i")

  points(seq(365, 1100, 15), corM[,1], t = "l", col = col.main[1], lwd = 2)
  polygon(c(seq(365, 1100, 15), rev(seq(365, 1100, 15))), c(confL[,1], rev(confU[,1])),col = col.alpha2[1], border=NA)

  points(seq(365, 1100, 15), corM[,2], t = "l", col = col.main[2], lwd = 2)
  polygon(c(seq(365, 1100, 15), rev(seq(365, 1100, 15))), c(confL[,2], rev(confU[,2])),col = col.alpha2[2], border=NA)

  axis(2, at = seq(-.5, .1, .1), las=2)
  axis(2, at = seq(-.5, .1, .05), rep("", length(seq(-.5, .1, .05))), las=2)
  axis(1, at = seq(100, 1100, 100), rep("",length(seq(100, 1100, 100))))

  abline(h = -0.5 + (2.75/50), lty=1, col = "lightgray")
  points(seq(365, 1100, 15), (lods / 50) + -0.5, t = "l", lwd=2)
  axis(4, at = c(-0.5 + 0.04, -0.5 + 0.08, -0.5 + 0.12, -0.5 + 0.16), c(2,4,6,8), las = 2,  cex.axis =1.0)


  legend("topleft", c("Female", "Male"), col = col.main, lwd=4, bg = "white", ncol=4, bty = "n")
dev.off()

### Continue


bCor <- cor(cdata[, "adjLongevity"], cdata[, "adjBw12"], use = "pair", method = "spearman")
bCor.f <- cor(cdata[which(cdata[, "sex"] == 0), "adjLongevity"], cdata[which(cdata[, "sex"] == 0), "adjBw12"], use = "pair", method = "spearman")
bCor.m <- cor(cdata[which(cdata[, "sex"] == 1), "adjLongevity"], cdata[which(cdata[, "sex"] == 1), "adjBw12"], use = "pair", method = "spearman")

computeDiffCor <- function(mcross, gtsp, cdata, sex = c(0, 1), method = "pearson"){
  allCor <- c()
  allN <- c()
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsp[, grep(marker, colnames(gtsp))]
    gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
      if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
    }))
    CH <- which(gts == "AC" & cdata[, "sex"] %in% sex)
    BH <- which(gts == "BC" & cdata[, "sex"] %in% sex)
    CD <- which(gts == "AD" & cdata[, "sex"] %in% sex)
    BD <- which(gts == "BD" & cdata[, "sex"] %in% sex)
    cCH <- NA; cBH <- NA; cCD <- NA; cBD <- NA
    if(length(CH) > 100) cCH <- cor(cdata[CH, "adjLongevity"], cdata[CH, "adjBw12"], use = "pair", method = method);
    if(length(BH) > 100) cBH <- cor(cdata[BH, "adjLongevity"], cdata[BH, "adjBw12"], use = "pair", method = method);
    if(length(CD) > 100) cCD <- cor(cdata[CD, "adjLongevity"], cdata[CD, "adjBw12"], use = "pair", method = method);
    if(length(BD) > 100) cBD <- cor(cdata[BD, "adjLongevity"], cdata[BD, "adjBw12"], use = "pair", method = method);

    allCor <- rbind(allCor, c(cCH, cBH, cCD, cBD))
    allN <- rbind(allN, c(length(CH), length(BH), length(CD), length(BD)))
  }
  rownames(allCor) <- colnames(pull.geno(mcross))
  colnames(allCor) <- c("CH", "BH", "CD", "BD")

  rownames(allN) <- colnames(pull.geno(mcross))
  colnames(allN) <- c("CH", "BH", "CD", "BD")

  # All correlations are NA, remove marker
  isN <- which(apply(apply(allCor,1,is.na), 2, sum) == 4)
  if (length(isN) > 0) {
    allCor <- allCor[-isN,]
    allN <- allN[-isN,]
  }

  # At least 500 individuals
  ii <- which(apply(allN, 1, sum) < 500)
  if (length(ii) > 0) {
    allCor <- allCor[-ii,]
    allN <- allN[-ii,]
  }
  return(list(allCor, allN))
}

res.c <- computeDiffCor(mcross, gtsp, cdata, c(0,1), method = "spearman")
res.f <- computeDiffCor(mcross, gtsp, cdata, 0, method = "spearman")
res.m <- computeDiffCor(mcross, gtsp, cdata, 1, method = "spearman")


toP <- function(allCor, allN, bCor){
  ## correlation differences to P-value / LOD scores
  pM <- c()
  pC <- c()
  for(x in 1:nrow(allCor)){
    pvs <- c()

    # Compute significance of allele cor relative to the base correlation
    for(y in 1:ncol(allCor)){
      cor <- c(allCor[x, y], bCor)
      z <- .5*log((1.0 + cor)/(1.0 - cor))
      n <- c(allN[x, y], nrow(cdata))
      df <- n-3
      sumOfSq <- sum(df * z^2)
      sqOfSum <- sum(df * z)
      denom <- sum(df)
      Cv <- sumOfSq - (sqOfSum^2/ denom)
      pvs <- c(pvs, pchisq(Cv, 1, 0, FALSE))
    }
    pM <- rbind(pM, pvs)

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
  return(list(pC,pM))
}

p.c <- toP(res.c[[1]], res.c[[2]], bCor)
p.f <- toP(res.f[[1]], res.f[[2]], bCor.f)
p.m <- toP(res.m[[1]], res.m[[2]], bCor.m)


map <- cbind(Chr = unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1)), 
             Pos = as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2))))
rownames(map) <- colnames(pull.geno(mcross))

subset <- map
subset <- cbind(subset, cpos = NA)

gap <- 30000000
chr.start <- c(0)
chr.mids <- c()
cp <- 0
for(chr in c(1:19, "X")){
  cl <- max(as.numeric(subset[which(subset[,1] == chr),2]))
  chr.start <- c(chr.start, cl + cp + gap)
  subset[which(subset[,1] == chr), "cpos"] <- as.numeric(subset[which(subset[,1] == chr), 2]) + cp
  if(chr == "X"){ chr <- 20 }
  cat(chr, " ", cl, " ", cp, " ", gap, " ", chr.start[chr], "\n")
  chr.mids <- c(chr.mids, chr.start[as.numeric(chr)] + cl/2)
  cp = cl + cp + gap
}


setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/__Tables")
write.table(cbind(round(-log10(p.c[[1]]),2), round(res.c[[1]],2)), file = "CTL_BW12_T365_C.txt", sep="\t", quote = FALSE)
write.table(cbind(round(-log10(p.m[[1]]),2), round(res.m[[1]],2)), file = "CTL_BW12_T365_M.txt", sep="\t", quote = FALSE)
write.table(cbind(round(-log10(p.f[[1]]),2), round(res.f[[1]],2)), file = "CTL_BW12_T365_F.txt", sep="\t", quote = FALSE)

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")
pdf(paste0("CTL_mapping_12mo.pdf"), width = 36, height = 12)
par(cex=2)
par(cex.axis=1.5)
plot(c(0, max(chr.start)), y = c(0, 12), t = 'n', ylab = "LOD", xlab = "Chromosome",xaxt="n", las=2, main = "CTL: Longevity (> 365 days) x Bodyweight 12 Months")
for(x in c(1:19, "X")) {
  points(subset[which(subset[,1] == x),"cpos"], -log10(p.c[[1]][rownames(subset)[which(subset[,1] == x)]]), t = 'l', col = "black",lwd=2, pch=20)
  points(subset[which(subset[,1] == x),"cpos"], -log10(p.f[[1]][rownames(subset)[which(subset[,1] == x)]]), t = 'l', col = "#FF3333",lwd=2, pch=20)
  points(subset[which(subset[,1] == x),"cpos"], -log10(p.m[[1]][rownames(subset)[which(subset[,1] == x)]]), t = 'l', col = "#00AEEF",lwd=2, pch=20)
}
abline(h = 2.75, lty=2, col = "darkolivegreen2", lwd=2)
abline(h = 4, lty=2, col = "darkolivegreen4", lwd=2)
axis(1, at = chr.mids, paste0("", c(1:19, "X")), cex.axis=1.2, las=1)
legend("topleft", c("FDR 5%", "FDR 1%"), lwd=2, lty=c(2,2), col = c("darkolivegreen2", "darkolivegreen4"))
legend("topright", c("Combined", "Females", "Males"), lwd=2, col = c("black", "#FF3333", "#00AEEF"))
dev.off()


