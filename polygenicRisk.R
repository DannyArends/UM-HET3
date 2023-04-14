setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

regions <- read.table("regions_4way_sorted_final.txt", sep = "\t", header=TRUE)
regions[,"Top"] <- gsub(",", "", regions[,"Top"])

as.DNI <- function(regions, col = "BALB.C3H"){
  bm <- as.numeric(unlist(lapply(strsplit(regions[, col], " ± "), "[",1)))
  bs <- as.numeric(unlist(lapply(strsplit(regions[, col], " ± "), "[",2)))

  isN <- abs(bm) <= abs(bs) # NEUTRAL
  isI <- abs(bm) >= abs(bs) & bm > 0 #Increasing
  isD <- abs(bm) >= abs(bs) & bm < 0 #Decreasing
  isV <- as.numeric(isD) * 1 + as.numeric(isN) * 2 + as.numeric(isI) * 3
  return(c("D", "N", "I")[isV])
}

risks <- cbind(AC = as.DNI(regions, "BALB.C3H"), AD = as.DNI(regions, "BALB.DBA"), BC = as.DNI(regions, "B6J.C3H"), BD = as.DNI(regions, "B6J.DBA"))
risks <- cbind(regions[, "Chr"],regions[, "Top"], risks)

risks2 <- risks

risks2[, "AC"] <-  as.numeric(unlist(lapply(strsplit(regions[, "BALB.C3H"], " ± "), "[",1)))
risks2[, "AD"] <-  as.numeric(unlist(lapply(strsplit(regions[, "BALB.DBA"], " ± "), "[",1)))
risks2[, "BC"] <-  as.numeric(unlist(lapply(strsplit(regions[, "B6J.C3H"], " ± "), "[",1)))
risks2[, "BD"] <-  as.numeric(unlist(lapply(strsplit(regions[, "B6J.DBA"], " ± "), "[",1)))

risks2[risks == "N"] <- NA

# I = Increase, D = Decrease, N = Neutral
#risks <- read.table("risks.txt", sep = "\t", header=TRUE)
# Effect size
#risks2 <- read.table("risk2.txt", sep = "\t", header=TRUE)

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
  sum(as.numeric(x), na.rm=TRUE)
})

#plot(pull.pheno(mcross)[iix, "longevity"] ~ mRp2[iix])



layout(matrix(c(1,1,2,2,2,1,1,2,2,2), 2, 5, byrow = TRUE))

plot(c(1, 13), c(750, 1000), t = 'n', xaxt='n', xlab = "(Inc - Dec) alleles", ylab="Lifespan", las=2, yaxs='i')
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
points(aa$conf[1,], pch = "-", cex=2, t = "h", col = xx, lwd=4)
points(aa$conf[1,], pch = "-", cex=3)
axis(1, at = 1:13, seq(-6,6,1))


mRp2c <- cut(mRp2, breaks = seq(-150, 150, 25), labels = FALSE)
iix <- which(pull.pheno(mcross)[, "longevity"] > 365)
aa <- boxplot(pull.pheno(mcross)[iix, "longevity"] ~ mRp2c[iix], plot=FALSE)

xx <- colorRampPalette(c("#E41A1C", "white", "#377EB8"))( 12 )

plot(c(1, 12), c(750, 1000), t = 'n', xaxt='n', xlab = "Expected effect on longevity ", ylab="Lifespan", las=2, yaxs='i')

phe <- pull.pheno(mcross)[iix, "longevity"]
sds <- c()
for(x in seq(1,12,1)){ sds <- c(sds, sd(phe[which(mRp2c[iix] == x)]) / sqrt(length(phe[which(mRp2c[iix] == x)]))); }


rect(1:12 - 0.25, rep(750,18), 1:12 + 0.25, aa$stats[3,], col = xx)
#points(aa$stats[3,], pch = 19, t = 'b')
#points(aa$stats[3,] + sds, pch = "-", cex=2, t = "h", lwd=2)
#points(aa$stats[3,] + sds, pch = "-", cex=3)
#points(aa$stats[3,] - sds, pch = "-", cex=2, t = "h", col = xx, lwd=2)
#points(aa$stats[3,] - sds, pch = "-", cex=3)
points(aa$conf[2,], pch = "-", cex=2, t = "h", lwd=2)
points(aa$conf[2,], pch = "-", cex=3)
points(aa$conf[1,], pch = "-", cex=2, t = "h", col = xx, lwd=4)
points(aa$conf[1,], pch = "-", cex=3)
mRp2c <- cut(mRp2, breaks = seq(-150, 150, 10))

axis(1, at = 1:18, names(table(mRp2c))[7:(7+17)])


