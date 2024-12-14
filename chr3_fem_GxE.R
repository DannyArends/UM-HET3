setwd("C:/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Github/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)
nxo <- countXO(mcross)

cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                    bw6 = pull.pheno(mcross)[, "bw6"], 
                    bw12 = pull.pheno(mcross)[, "bw12"], 
                    bw18 = pull.pheno(mcross)[, "bw18"], 
                    bw24 = pull.pheno(mcross)[, "bw24"], 
                    nxo = as.numeric(nxo), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
idx <- which(cdata[, "longevity"] > 365 & cdata[, "sex"] == 0)

cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lods.c <- c()
lm.null <- lm(longevity ~ site + treatment + cohort + 0, data = cdata)
for(marker in colnames(pull.geno(mcross))){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  lm.alt <- lm(longevity ~ site + treatment + cohort + mp + 0, data = cdata)
  n <- sum(!is.na(lm.alt$resid))
  lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
  lods.c <- c(lods.c, lod)
}
names(lods.c) <- colnames(pull.geno(mcross))


# Plot the QTL profile
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))
chr.l <- c()
for(chr in unique(chrs)){ chr.l <- c(chr.l, max(positions[which(chrs==chr)])); }
names(chr.l) <- unique(chrs)


plot(c(-1000000, 5000000 + sum(chr.l) + 19 * 25000000), c(0,8), main = "QTL profiles for longevity", t = 'n', xaxt='n', xlab="Chromosome", ylab="LOD", xaxs="i", yaxs="i", las=2)
#abline(h = seq(0,10,2), lty=2)
s <- 0
h <- c()
for(chr in unique(chrs)){
  rect(s, 0, s + chr.l[chr] + 5000000, 12, col = rgb(0.9,0.9,0.9,0.5),border=NA)
  points(positions[chrs==chr] + s, lods.c[chrs == chr], t = 'l', lwd=3, col="black")
  h <- c(h, s + (chr.l[chr]/2) + 25000000)
  s <- s + chr.l[chr] + 25000000
}
axis(1, at = h, unique(chrs))
abline(h = c(4.25, 4.95), col = c("orange", "green"), lty=2)
legend("topright", c("Female controls + NDE"), lwd=c(3), col=c("black"), bg="white")
legend("topleft", c("p = 0.05", "p = 0.01"), lwd=1, lty=2, col=c("orange", "green"), bg="white")
box()

marker <- "3_103803615"
mp <- gtsp[, grep(marker, colnames(gtsp))]
mlm <- lm(longevity ~ site + treatment + cohort + mp + mp:treatment + 0, data = cdata)
anova(mlm)

vioplot(residuals(mlm) + mean(cdata[, "longevity"]) ~ apply(mp,1,which.min), xaxt='n', xlab = "GT", main = "CTRL Females at Chr3 top")
axis(1, at = 1:4, colnames(mp))
