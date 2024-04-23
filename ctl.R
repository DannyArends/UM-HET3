setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
mapAll <- scanone(mcross)[,1:2]

genotypes <- pull.geno(fill.geno(mcross, method = "maxmarginal", error.prob=0.01, min.prob = 0.9))

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    adjLongevity = NA, 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw6 = as.numeric(pull.pheno(mcross)[, "bw6"]),
                    adjBw6 = NA,
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
idx <- which(cdata[, "longevity"] > 365 & !is.na(cdata[, "bw6"]))
cdata <- cdata[idx,]
gtsM <- genotypes[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),0)

lm.null.bw6 <- lm(bw6 ~ sex + site + cohort, data = cdata)
cdata[, "adjBw6"] <- round(as.numeric(coef(lm.null.bw6)["(Intercept)"]) + residuals(lm.null.bw6),0)

bCor <- cor(cdata[, "adjLongevity"], cdata[, "adjBw6"], use = "pair", method = "spearman")

allCor <- c()
allN <- c()
for(x in 1:ncol(genotypes)){
  balb <- which(genotypes[,x] == 1)
  b6 <- which(genotypes[,x] == 2)
  c3h <- which(genotypes[,x] == 3)
  dba <- which(genotypes[,x] == 4)
  cBalb <- NA; cB6 <- NA; cC3h <- NA; cDba <- NA
  if(length(balb) > 100) cBalb <- cor(cdata[balb, "adjLongevity"], cdata[balb, "adjBw6"], use = "pair", method = "spearman");
  if(length(b6) > 100) cB6 <- cor(cdata[b6, "adjLongevity"], cdata[b6, "adjBw6"], use = "pair", method = "spearman");
  if(length(c3h) > 100) cC3h <- cor(cdata[c3h, "adjLongevity"], cdata[c3h, "adjBw6"], use = "pair", method = "spearman");
  if(length(dba) > 100) cDba <- cor(cdata[dba, "adjLongevity"], cdata[dba, "adjBw6"], use = "pair", method = "spearman");
  allCor <- rbind(allCor, c(cBalb, cB6, cC3h, cDba))
  allN <- rbind(allN, c(length(balb), length(b6), length(c3h), length(dba)))
}
rownames(allCor) <- colnames(genotypes)
colnames(allCor) <- c("BALB/cByJ | C3H/HeJ", "C57BL/6J | C3H/HeJ", "BALB/cByJ | DBA/2J", "C57BL/6J | DBA/2J")

rownames(allN) <- colnames(genotypes)
colnames(allN) <- c("BALB/cByJ | C3H/HeJ", "C57BL/6J | C3H/HeJ", "BALB/cByJ | DBA/2J", "C57BL/6J | DBA/2J")


allCor <- allCor[-which(apply(apply(allCor,1,is.na),2,sum) == 4),]
write.table(allCor, file="correlations_Long_BW6.txt",sep="\t", quote=FALSE)

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
colz.c <- colorRampPalette(c("#41ab5d", "#fd8d3c", "red"))(30)

mids <- c(0, getHalf(as.numeric(pos[iix])), max(as.numeric(pos[iix])))

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
pdf(paste0("Figure_3_CTL_Chr6.pdf"), width = 16, height = 8)

plot(c(0, max(as.numeric(pos[iix]))), y = c(1,5), t = "n", yaxs= "i", xaxs= "i",xlab = "Position (Mb)", ylab = "Genotype", xaxt="n", yaxt = "n", main = "CTL Chromosome 6")
i <- 1
for(x in 1:nrow(allCor[iix,])){
  rect(mids[i], 1, mids[i+1], 2, col = colz.c[round(abs(allCor[iix[x],1]) * 100)], border=colz.c[round(abs(allCor[iix[x],1]) * 100)])
  rect(mids[i], 2, mids[i+1], 3, col = colz.c[round(abs(allCor[iix[x],3]) * 100)], border=colz.c[round(abs(allCor[iix[x],3]) * 100)])
  rect(mids[i], 3, mids[i+1], 4, col = colz.c[round(abs(allCor[iix[x],2]) * 100)], border=colz.c[round(abs(allCor[iix[x],2]) * 100)])
  rect(mids[i], 4, mids[i+1], 5, col = colz.c[round(abs(allCor[iix[x],4]) * 100)], border=colz.c[round(abs(allCor[iix[x],4]) * 100)])
  i <- i + 1
}
axis(1, at = seq(0, max(as.numeric(pos[iix])), 10000000), seq(0, max(as.numeric(pos[iix])), 10000000) / 1000000)
axis(2, at = 0.5 + 1:4, c("CH", "CD", "BH", "BD"), las=2)
box()
dev.off()

## Chromosome 2, 6, and 10
pvs <- c()
for(x in 1:nrow(allCor)){
  cor <- allCor[x, ] # ["2_24606612",]
  z <- .5*log((1.0 + cor)/(1.0 - cor))
  n <- allN[x, ] #["2_24606612",]
  df <- n-3
  sumOfSq <- sum(df * z^2)
  sqOfSum <- sum(df * z)
  denom <- sum(df)
  Cv <- sumOfSq - (sqOfSum^2/ denom)
  pvs <- c(pvs, pchisq(Cv, 1, 0, FALSE))
}

names(pvs) <- rownames(allCor)

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

subset <- map[which(map[,1] %in% c(1:19, "X")),]
subset <- subset[which(rownames(subset) %in% rownames(allCor)),]
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

plot(c(0, max(chr.start)), y = c(-0.3, 0), t = 'n', ylab = "COR", xlab = "Chromosome",xaxt="n", las=2, main = "CTL Longevity x BW6")
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

plot(-log10(pvs), t = 'l')


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


