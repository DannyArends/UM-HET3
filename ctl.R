setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

#
# Use a different clor scheme (blue - red)
#

### TODO: CTL for Males and Females separately

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

bCor <- cor(cdata[, "adjLongevity"], cdata[, "adjBw6"], use = "pair")
bCor.f <- cor(cdata[which(cdata[, "sex"] == 0), "adjLongevity"], cdata[which(cdata[, "sex"] == 0), "adjBw6"], use = "pair")
bCor.m <- cor(cdata[which(cdata[, "sex"] == 1), "adjLongevity"], cdata[which(cdata[, "sex"] == 1), "adjBw6"], use = "pair")


allCor <- c()
allN <- c()
for(marker in colnames(pull.geno(mcross))){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  CH <- which(gts == "AC")
  BH <- which(gts == "BC")
  CD <- which(gts == "AD")
  BD <- which(gts == "BD")
  cCH <- NA; cBH <- NA; cCD <- NA; cBD <- NA
  if(length(CH) > 100) cCH <- cor(cdata[CH, "adjLongevity"], cdata[CH, "adjBw6"], use = "pair");
  if(length(BH) > 100) cBH <- cor(cdata[BH, "adjLongevity"], cdata[BH, "adjBw6"], use = "pair");
  if(length(CD) > 100) cCD <- cor(cdata[CD, "adjLongevity"], cdata[CD, "adjBw6"], use = "pair");
  if(length(BD) > 100) cBD <- cor(cdata[BD, "adjLongevity"], cdata[BD, "adjBw6"], use = "pair");

  allCor <- rbind(allCor, c(cCH, cBH, cCD, cBD))
  allN <- rbind(allN, c(length(CH), length(BH), length(CD), length(BD)))
}
rownames(allCor) <- colnames(pull.geno(mcross))
colnames(allCor) <- c("CH", "BH", "CD", "BD")

rownames(allN) <- colnames(pull.geno(mcross))
colnames(allN) <- c("CH", "BH", "CD", "BD")

isN <- which(apply(apply(allCor,1,is.na),2,sum) == 4)

allCor <- allCor[-isN,]
allN <- allN[-isN,]

ii <- which(apply(allN,1,sum) < 1000)

allCor <- allCor[-ii,]
allN <- allN[-ii,]

write.table(allCor, file="correlations_Long_BW6.txt",sep="\t", quote=FALSE)

## correlation differences to P-value / LOD scores
pM <- c()
pC <- c()
for(x in 1:nrow(allCor)){
  pvs <- c()
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

  cor <- allCor[x,]
  z <- .5*log((1.0 + cor)/(1.0 - cor))
  n <- allN[x,]
  df <- n-3
  sumOfSq <- sum(df * z^2)
  sqOfSum <- sum(df * z)
  denom <- sum(df)
  Cv <- sumOfSq - (sqOfSum^2/ denom)
  pC <- c(pC, pchisq(Cv, 1, 0, FALSE))
}
names(pC) <- rownames(allCor)

map <- cbind(Chr = unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1)), 
             Pos = as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2))))
rownames(map) <- colnames(pull.geno(mcross))

subset <- map[which(rownames(map) %in% rownames(allCor)),]
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

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
pdf(paste0("Figure_3_CTL_mapping.pdf"), width = 36, height = 12)
par(cex=2)
par(cex.axis=1.5)
plot(c(0, max(chr.start)), y = c(0, 8), t = 'n', ylab = "LOD", xlab = "Chromosome",xaxt="n", las=2, main = "CTL: Longevity (≥ 365 days) x Bodyweight 6 Months")
for(x in c(1:19, "X")){
  points(subset[which(subset[,1] == x),"cpos"], -log10(pC[rownames(subset)[which(subset[,1] == x)]]), t = 'l', col = "black",lwd=2)
}
abline(h = 2.75, lty=2, col = "darkolivegreen2", lwd=2)
abline(h = 4, lty=2, col = "darkolivegreen4", lwd=2)
axis(1, at = chr.mids, paste0("", c(1:19, "X")), cex.axis=1.2, las=1)
legend("topleft", c("FDR 5%", "FDR 1%"), lwd=2, lty=c(2,2), col = c("darkolivegreen2", "darkolivegreen4"))
dev.off()

write.table(round(allCor[names(which(p.adjust(pC, "fdr") < 0.05)),],2), "CTL_results.txt", sep = "\t")


round(allCor[names(which(p.adjust(pC, "fdr") < 0.05)),],2)

round(-log10(pC),2)

all <- c("1_3010274", "2_60201233", "3_156416491", "4_11234654", "7_16072018", "8_95039510", "13_53167285", "X_105166309")
names(all) <- c("CTL1a", "CTL2a", "CTL3a", "CTL4a", "CTL7a", "CTL8a", "CTL13a", "CTLXa")

genotypes <- c()
for(marker in all){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  genotypes <- cbind(genotypes, gts)
}
colnames(genotypes) <- all


for(i in 1:length(all)){
  x <- all[i]
  CH <- which(genotypes[,x] == "AC")
  BH <- which(genotypes[,x] == "BC")
  CD <- which(genotypes[,x] == "AD")
  BD <- which(genotypes[,x] == "BD")
  
  col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
  col.main <- adjustcolor( col.main, alpha.f = 0.6)

  pdf(paste0("Figure_3_CTL_",names(all)[i],".pdf"), width = 8, height = 8)
  plot(c(20, 55), c(365, 1250), t = "n", xlab = "Adjusted bodyweight (6mo)", ylab = "Adjusted longevity", main = paste0("CTL @ ", x), log = "y", yaxt = "n")
  points(cdata[CH, "adjBw6"], cdata[CH, "adjLongevity"], pch = 19, col = col.main[1])
  points(cdata[CD, "adjBw6"], cdata[CD, "adjLongevity"], pch = 19, col = col.main[2])
  points(cdata[BH, "adjBw6"], cdata[BH, "adjLongevity"], pch = 19, col = col.main[3])
  points(cdata[BD, "adjBw6"], cdata[BD, "adjLongevity"], pch = 19, col = col.main[4])

  abline(lm(cdata[CH, "adjLongevity"] ~ cdata[CH, "adjBw6"]), col = col.main[1], lwd=2, untf=T)
  abline(lm(cdata[CD, "adjLongevity"] ~ cdata[CD, "adjBw6"]), col = col.main[2], lwd=2, untf=T)
  abline(lm(cdata[BH, "adjLongevity"] ~ cdata[BH, "adjBw6"]), col = col.main[3], lwd=2, untf=T)
  abline(lm(cdata[BD, "adjLongevity"] ~ cdata[BD, "adjBw6"]), col = col.main[4], lwd=2, untf=T)

  axis(2, at = seq(400, 1200,200), seq(400, 1200,200), las=2)

  legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=19, bg = "white", ncol=4, bty = "n")
  dev.off()
}


#### OLD code there be Dragons down here

pvs[names(apply(!apply(genotypes[,names(which(apply(allCor,1, function(x){any(x > -0.1)})))],2,is.na),2,sum))]

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

subset <- map[which(map[,1] %in% c(1:19, "X")),]
subset <- subset[which(rownames(subset) %in% rownames(allCor)),]
subset <- cbind(subset, cpos = NA)
gap <- 10000000
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




chr <- unlist(lapply(strsplit(rownames(allCor), "_"), "[",1))
pos <- unlist(lapply(strsplit(rownames(allCor), "_"), "[",2))

getHalf <- function(x){
  mids <- c()
  cur <- x[1]
  for(y in 2:length(x)){
    mids <- c(mids, (cur + x[y]) / 2)
    cur <- x[y]
  }
  return(mids)
}

library(RColorBrewer)
library(svglite)
iix <- which(chr == 6)
pos[iix]
allCor[iix, ]
colz.c <- colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(30)

mids <- c(0, getHalf(as.numeric(pos[iix])), max(as.numeric(pos[iix])))

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
#pdf(paste0("Figure_3_CTL_Chr6.pdf"), width = 16, height = 8)

plot(c(0, max(as.numeric(pos[iix]))), y = c(1,5), t = "n", yaxs= "i", xaxs= "i",xlab = "Position (Mb)", ylab = "Genotype", xaxt="n", yaxt = "n", main = "CTL Chromosome 6")
i <- 1
for(x in 1:nrow(allCor[iix,])){
  rect(mids[i], 1, mids[i+1], 2, col = colz.c[round(abs(allCor[iix[x],1]) * 100)], border=colz.c[round(abs(allCor[iix[x],"CH"]) * 100)])
  rect(mids[i], 2, mids[i+1], 3, col = colz.c[round(abs(allCor[iix[x],3]) * 100)], border=colz.c[round(abs(allCor[iix[x],"CD"]) * 100)])
  rect(mids[i], 3, mids[i+1], 4, col = colz.c[round(abs(allCor[iix[x],2]) * 100)], border=colz.c[round(abs(allCor[iix[x],"BH"]) * 100)])
  rect(mids[i], 4, mids[i+1], 5, col = colz.c[round(abs(allCor[iix[x],4]) * 100)], border=colz.c[round(abs(allCor[iix[x],"BD"]) * 100)])
  i <- i + 1
}
axis(1, at = seq(0, max(as.numeric(pos[iix])), 10000000), seq(0, max(as.numeric(pos[iix])), 10000000) / 1000000)
axis(2, at = 0.5 + 1:4, c("CH", "CD", "BH", "BD"), las=2)
box()
#dev.off()


### Plotting correlation at the top markers
apply(!apply(genotypes[,names(which(apply(allCor,1, function(x){any(x > -0.1)})))],2,is.na),2,sum)

m2 <- "2_25363944"
m6 <- "6_51450974"
m10 <- "10_127745179"

for(x in c(m2, m6, m10)){
  CH <- which(genotypes[,x] == 1)
  BH <- which(genotypes[,x] == 2)
  CD <- which(genotypes[,x] == 3)
  BD <- which(genotypes[,x] == 4)
  
  col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
  col.main <- adjustcolor( col.main, alpha.f = 0.6)

  pdf(paste0("Figure_3_CTL_",x,".pdf"), width = 8, height = 8)
  plot(c(20, 55), c(365, 1250), t = "n", xlab = "Adjusted bodyweight (6mo)", ylab = "Adjusted longevity", main = paste0("CTL @ ", x), log = "y", yaxt = "n")
  points(cdata[CH, "adjBw6"], cdata[CH, "adjLongevity"], pch = 19, col = col.main[1])
  points(cdata[CD, "adjBw6"], cdata[CD, "adjLongevity"], pch = 19, col = col.main[2])
  points(cdata[BH, "adjBw6"], cdata[BH, "adjLongevity"], pch = 19, col = col.main[3])
  points(cdata[BD, "adjBw6"], cdata[BD, "adjLongevity"], pch = 19, col = col.main[4])

  abline(lm(cdata[CH, "adjLongevity"] ~ cdata[CH, "adjBw6"]), col = col.main[1], lwd=2, untf=T)
  abline(lm(cdata[CD, "adjLongevity"] ~ cdata[CD, "adjBw6"]), col = col.main[2], lwd=2, untf=T)
  abline(lm(cdata[BH, "adjLongevity"] ~ cdata[BH, "adjBw6"]), col = col.main[3], lwd=2, untf=T)
  abline(lm(cdata[BD, "adjLongevity"] ~ cdata[BD, "adjBw6"]), col = col.main[4], lwd=2, untf=T)

  axis(2, at = seq(400, 1200,200), seq(400, 1200,200), las=2)

  legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=19, bg = "white", ncol=4, bty = "n")
  dev.off()
}


### Mapping for P-values
pdf("Figure_3_CTLeffects.pdf", width = 36, height = 12)

plot(c(0, max(chr.start)), y = c(-0.3, 0), t = 'n', ylab = "Correlation", xlab = "Chromosome",xaxt="n", las=2, main = "CTL Longevity x BW6")
col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
for(x in c(1:19, "X")){
  onC <- which(subset[,1] == x)
  #points(subset[onC,"cpos"], -log10(pvs[rownames(subset)[onC]]), t = 'l', col = "black",lwd=2)
  points(subset[onC,"cpos"], allCor[rownames(subset)[onC],1], t = 'l', col = col.main[1],lwd=2)
  points(subset[onC,"cpos"], allCor[rownames(subset)[onC],2], t = 'l', col = col.main[2],lwd=2)
  points(subset[onC,"cpos"], allCor[rownames(subset)[onC],3], t = 'l', col = col.main[3],lwd=2)
  points(subset[onC,"cpos"], allCor[rownames(subset)[onC],4], t = 'l', col = col.main[4],lwd=2)

}
abline(h = 3.65, lty=2, col = "green")
axis(1, at = chr.mids, paste0("", c(1:19, "X")), cex.axis=1.2, las=1)
legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=19, bg = "white", ncol=4, bty = "n")
dev.off()


plot(-log10(pvs), t = 'l')



### Old figure

library(RColorBrewer)
colz <- colorRampPalette(c("white", "red"))(15)[15:1]
op <- par(mar = c(4,10,3,1))
layout(t(matrix(c(1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,0), 11, 2, byrow = FALSE)))
layout.show(2)

image(1:nrow(allCor), 1:4, allCor, breaks = seq(0.0, -0.3, -0.02), col = colz, xaxt='n', yaxt='n', ylab="", xlab="Marker", main = "Correlation between longevity and body weight at 6 month within genotype")
nMchr <- table(unlist(lapply(strsplit(rownames(allCor), "_"),"[",1)))
nMchr <- nMchr[c(1:19,"X")]
cPos <- 1
mPoss <- c()
for(chr in 1:length(nMchr)){
  mPoss <- c(mPoss, cPos + nMchr[chr] / 2)
  cPos <- cPos + nMchr[chr]
  abline(v = cPos +0.5, col = "white", lwd=2)
}
axis(1, at = mPoss, c(1:19, "X"))
axis(2, at = 1:4, colnames(allCor), las=2)
abline(h = 1:4 +0.5, col = "white", lwd=2)
box()

image(t(rev(seq(0.0, -0.3, -0.02))), breaks = seq(0.0, -0.3, -0.02), col = colz, xaxt='n', yaxt='n', ylab="")
axis(2, at = seq(0,1,1/6), rev(seq(0.0, -0.3, -0.05)), las=2)

x <- "2_25460620"
g1 <- which(genotypes[,x] == 1)
g2 <- which(genotypes[,x] == 2)
g3 <- which(genotypes[,x] == 3)
g4 <- which(genotypes[,x] == 4)

plot(cdata[, "adjLongevity"], cdata[, "adjBw6"], col = genotypes[,x], pch=19)

write.table(allCor[names(which(apply(apply(allCor,1,function(x){x > -0.1}),2,any))),], file="corD.txt",sep="\t", quote=FALSE)


