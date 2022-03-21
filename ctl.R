setwd("D:/Ddrive/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Github/UM-HET3/files")

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
}
rownames(allCor) <- colnames(genotypes)
colnames(allCor) <- c("BALB/cByJ | C3H/HeJ", "C57BL/6J | C3H/HeJ", "BALB/cByJ | DBA/2J", "C57BL/6J | DBA/2J")

allCor <- allCor[-which(apply(apply(allCor,1,is.na),2,sum) == 4),]

library(RColorBrewer)
colz <- colorRampPalette(c("blue", "red"))(15)[15:1]
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


