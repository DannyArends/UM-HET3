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
idx <- which(cdata[, "longevity"] >= 0 & !is.na(cdata[, "bw6"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw6 <- lm(bw6 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw6"] <- round(as.numeric(coef(lm.null.bw6)["(Intercept)"]) + residuals(lm.null.bw6), 2)

bCor <- cor(cdata[, "adjLongevity"], cdata[, "adjBw6"], use = "pair", method = "spearman")
bCor.f <- cor(cdata[which(cdata[, "sex"] == 0), "adjLongevity"], cdata[which(cdata[, "sex"] == 0), "adjBw6"], use = "pair", method = "spearman")
bCor.m <- cor(cdata[which(cdata[, "sex"] == 1), "adjLongevity"], cdata[which(cdata[, "sex"] == 1), "adjBw6"], use = "pair", method = "spearman")

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
    if(length(CH) > 100) cCH <- cor(cdata[CH, "adjLongevity"], cdata[CH, "adjBw6"], use = "pair", method = method);
    if(length(BH) > 100) cBH <- cor(cdata[BH, "adjLongevity"], cdata[BH, "adjBw6"], use = "pair", method = method);
    if(length(CD) > 100) cCD <- cor(cdata[CD, "adjLongevity"], cdata[CD, "adjBw6"], use = "pair", method = method);
    if(length(BD) > 100) cBD <- cor(cdata[BD, "adjLongevity"], cdata[BD, "adjBw6"], use = "pair", method = method);

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
    notNA <- which(!is.na(allCor[x,]))
    cor <- allCor[x, notNA]
    z <- .5*log((1.0 + cor)/(1.0 - cor))
    n <- allN[x,notNA]
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


setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")
pdf(paste0("CTL_mapping_42day.pdf"), width = 36, height = 12)
par(cex=2)
par(cex.axis=1.5)
plot(c(0, max(chr.start)), y = c(0, 8), t = 'n', ylab = "LOD", xlab = "Chromosome",xaxt="n", las=2, main = "CTL: Longevity (≥ 42 days) x Bodyweight 42 Days")
for(x in c(1:19, "X")) {
  points(subset[which(subset[,1] == x),"cpos"], -log10(p.c[[1]][rownames(subset)[which(subset[,1] == x)]]), t = 'b', col = "black",lwd=2, pch=20)
  points(subset[which(subset[,1] == x),"cpos"], -log10(p.f[[1]][rownames(subset)[which(subset[,1] == x)]]), t = 'b', col = "hotpink",lwd=2, pch=20)
  points(subset[which(subset[,1] == x),"cpos"], -log10(p.m[[1]][rownames(subset)[which(subset[,1] == x)]]), t = 'b', col = "blue",lwd=2, pch=20)
}
abline(h = 2.75, lty=2, col = "darkolivegreen2", lwd=2)
abline(h = 4, lty=2, col = "darkolivegreen4", lwd=2)
axis(1, at = chr.mids, paste0("", c(1:19, "X")), cex.axis=1.2, las=1)
legend("topleft", c("FDR 5%", "FDR 1%"), lwd=2, lty=c(2,2), col = c("darkolivegreen2", "darkolivegreen4"))
legend("topright", c("Combined", "Females", "Males"), lwd=2, col = c("black", "hotpink", "blue"))
dev.off()

round(-log10(p.f[[1]],2)



