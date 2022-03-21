setwd("D:/Ddrive/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Github/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)


# I = Increase, D = Decrease, N = Neutral
risks <- read.table("risks.txt", sep = "\t", header=TRUE)
# Effect size
risks2 <- read.table("risk2.txt", sep = "\t", header=TRUE)

mRm <- matrix(NA, nrow(gtsp), nrow(risks),dimnames = list(rownames(gtsp),  paste0(risks[,1], "_", risks[,2])))
mRm2 <- matrix(NA, nrow(gtsp), nrow(risks),dimnames = list(rownames(gtsp),  paste0(risks[,1], "_", risks[,2])))

for(x in 1:nrow(risks)){
  marker <- paste0(risks[x,1], "_", risks[x,2])
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  mRm[which(gts == "AC"),x] <- risks[x, "AC"]
  mRm[which(gts == "AD"),x] <- risks[x, "AD"]
  mRm[which(gts == "BC"),x] <- risks[x, "BC"]
  mRm[which(gts == "BD"),x] <- risks[x, "BD"]
  
  mRm2[which(gts == "AC"),x] <- risks2[x, "AC"]
  mRm2[which(gts == "AD"),x] <- risks2[x, "AD"]
  mRm2[which(gts == "BC"),x] <- risks2[x, "BC"]
  mRm2[which(gts == "BD"),x] <- risks2[x, "BD"]
}

mRp <- apply(mRm,1, function(x){
  length(which(x == "I")) - length(which(x == "D"))
})

mRp2 <- apply(mRm2,1, function(x){
  sum(x, na.rm=TRUE)
})

plot(pull.pheno(mcross)[iix, "longevity"] ~ mRp2[iix])



layout(matrix(c(1,1,2,2,2,1,1,2,2,2), 2, 5, byrow = TRUE))

plot(c(1, 13), c(750, 1000), t = 'n', xaxt='n', xlab = "(Inc - Dec) alleles", ylab="Longevity", las=2, yaxs='i')
iix <- which(mRp > -7 & mRp < 7 & pull.pheno(mcross)[, "longevity"] > 365)
xx <- colorRampPalette(c("#E41A1C", "white", "#377EB8"))( 13 )
aa <- boxplot(pull.pheno(mcross)[iix, "longevity"] ~ mRp[iix], plot=FALSE)

phe <- pull.pheno(mcross)[iix, "longevity"]
sds <- c()
for(x in seq(-6,6,1)){ sds <- c(sds, sd(phe[which(mRp[iix] == x)]) / sqrt(length(phe[which(mRp[iix] == x)]))); }

rect(1:13 - 0.25, rep(750,13), 1:13 + 0.25, aa$stats[3,], col = xx)
#points(aa$stats[3,] + sds, pch = "-", cex=2, t = "h", lwd=2)
#points(aa$stats[3,] + sds, pch = "-", cex=3)
#points(aa$stats[3,] - sds, pch = "-", cex=2, t = "h", col = xx, lwd=2)
#points(aa$stats[3,] - sds, pch = "-", cex=3)
points(aa$conf[2,], pch = "-", cex=2, t = "h", lwd=2)
points(aa$conf[2,], pch = "-", cex=3)
points(aa$conf[1,], pch = "-", cex=2, t = "h", col = xx, lwd=2)
points(aa$conf[1,], pch = "-", cex=3)
axis(1, at = 1:13, seq(-6,6,1))


mRp2c <- cut(mRp2, breaks = seq(-150, 150, 10), labels = FALSE)
iix <- which(mRp2c > 6 & mRp2c < 25 & pull.pheno(mcross)[, "longevity"] > 365)
aa <- boxplot(pull.pheno(mcross)[iix, "longevity"] ~ mRp2c[iix], plot=FALSE)

xx <- colorRampPalette(c("#E41A1C", "white", "#377EB8"))( 18 )

plot(c(1, 18), c(750, 1000), t = 'n', xaxt='n', xlab = "Expected effect on longevity ", ylab="Longevity", las=2, yaxs='i')

phe <- pull.pheno(mcross)[iix, "longevity"]
sds <- c()
for(x in seq(7,24,1)){ sds <- c(sds, sd(phe[which(mRp2c[iix] == x)]) / sqrt(length(phe[which(mRp2c[iix] == x)]))); }


rect(1:18 - 0.25, rep(750,18), 1:18 + 0.25, aa$stats[3,], col = xx)
#points(aa$stats[3,], pch = 19, t = 'b')
#points(aa$stats[3,] + sds, pch = "-", cex=2, t = "h", lwd=2)
#points(aa$stats[3,] + sds, pch = "-", cex=3)
#points(aa$stats[3,] - sds, pch = "-", cex=2, t = "h", col = xx, lwd=2)
#points(aa$stats[3,] - sds, pch = "-", cex=3)
points(aa$conf[2,], pch = "-", cex=2, t = "h", lwd=2)
points(aa$conf[2,], pch = "-", cex=3)
points(aa$conf[1,], pch = "-", cex=2, t = "h", col = xx, lwd=2)
points(aa$conf[1,], pch = "-", cex=3)
mRp2c <- cut(mRp2, breaks = seq(-150, 150, 10))

axis(1, at = 1:18, names(table(mRp2c))[7:(7+17)])
