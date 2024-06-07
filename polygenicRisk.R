library(svglite)

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)
sex <- pull.pheno(mcross)[, "sex"]

# Include animals older than 365 days
above <- which(pull.pheno(mcross)[, "longevity"] > 365)

gtsp <- gtsp[above,]
sex <- sex[above]
observed <- pull.pheno(mcross)[above, "longevity"]
sex <- pull.pheno(mcross)[above, "sex"]

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3")

# Read the effects at top marker
regions <- read.table("regions_4way_merged_May24_effects.txt", sep="\t", header=TRUE, row.names=1)
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

# Risk as a locus count
risks <- cbind(AC = as.DNI(regions, "CH"), AD = as.DNI(regions, "CD"), BC = as.DNI(regions, "BH"), BD = as.DNI(regions, "BD"))
risks <- cbind(Chr = regions[, "Chr"], Top = regions[, "Top"], Affects = regions[, "Affects"], risks)

risks2 <- risks
# Risk as a count of effects
risks2[, "AC"] <-  as.numeric(unlist(lapply(strsplit(regions[, "CH"], " ± "), "[",1)))
risks2[, "AD"] <-  as.numeric(unlist(lapply(strsplit(regions[, "CD"], " ± "), "[",1)))
risks2[, "BC"] <-  as.numeric(unlist(lapply(strsplit(regions[, "BH"], " ± "), "[",1)))
risks2[, "BD"] <-  as.numeric(unlist(lapply(strsplit(regions[, "BD"], " ± "), "[",1)))
risks2[risks == "N"] <- NA

mRm <- matrix(NA, nrow(gtsp), nrow(risks), dimnames = list(rownames(gtsp),  paste0(risks[,1], "_", risks[,2])))
mRm2 <- matrix(NA, nrow(gtsp), nrow(risks), dimnames = list(rownames(gtsp),  paste0(risks[,1], "_", risks[,2])))

# Go through the risks/effects of each allele and sum
for(x in 1:nrow(risks)){
  marker <- paste0(risks[x,1], "_", risks[x,2])
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  if(risks[x, "Affects"] == "B"){ # Locus affects risk for both
    mRm[which(gts == "AC"),x] <- risks[x, "AC"]
    mRm[which(gts == "AD"),x] <- risks[x, "AD"]
    mRm[which(gts == "BC"),x] <- risks[x, "BC"]
    mRm[which(gts == "BD"),x] <- risks[x, "BD"]
  
    mRm2[which(gts == "AC"),x] <- risks2[x, "AC"]
    mRm2[which(gts == "AD"),x] <- risks2[x, "AD"]
    mRm2[which(gts == "BC"),x] <- risks2[x, "BC"]
    mRm2[which(gts == "BD"),x] <- risks2[x, "BD"]
  }

  if(risks[x, "Affects"] == "F"){ # Locus affects risk for females only
    mRm[which(gts == "AC" & sex == 0),x] <- risks[x, "AC"]
    mRm[which(gts == "AD" & sex == 0),x] <- risks[x, "AD"]
    mRm[which(gts == "BC" & sex == 0),x] <- risks[x, "BC"]
    mRm[which(gts == "BD" & sex == 0),x] <- risks[x, "BD"]
  
    mRm2[which(gts == "AC" & sex == 0),x] <- risks2[x, "AC"]
    mRm2[which(gts == "AD" & sex == 0),x] <- risks2[x, "AD"]
    mRm2[which(gts == "BC" & sex == 0),x] <- risks2[x, "BC"]
    mRm2[which(gts == "BD" & sex == 0),x] <- risks2[x, "BD"]
  }
  if(risks[x, "Affects"] == "M"){ # Locus affects risk for males only
    mRm[which(gts == "AC" & sex == 1),x] <- risks[x, "AC"]
    mRm[which(gts == "AD" & sex == 1),x] <- risks[x, "AD"]
    mRm[which(gts == "BC" & sex == 1),x] <- risks[x, "BC"]
    mRm[which(gts == "BD" & sex == 1),x] <- risks[x, "BD"]
  
    mRm2[which(gts == "AC" & sex == 1),x] <- risks2[x, "AC"]
    mRm2[which(gts == "AD" & sex == 1),x] <- risks2[x, "AD"]
    mRm2[which(gts == "BC" & sex == 1),x] <- risks2[x, "BC"]
    mRm2[which(gts == "BD" & sex == 1),x] <- risks2[x, "BD"]
  }
}

# Risk as increasing and decreasing alleles
mRp <- apply(mRm,1, function(x) { return(length(which(x == "I")) - length(which(x == "D"))) })

# Risk as sum off locus effects
mRp2 <- apply(mRm2,1, function(x){ return(sum(as.numeric(x), na.rm=TRUE)) })

#TODO (color by lifespan)

group <- c()
group.n <- c()
group.ac <- c()
for(x in sort(unique(mRp))){
  group.n <- c(group.n, length(which(mRp == x)))
  group <- c(group, round(mean(observed[which(mRp == x)]), 0))
  group.ac <- c(group.ac, x)
}

group <-group[which(group.n > 25)]
group.ac <-group.ac[which(group.n > 25)]

mRa <- mRp[which(mRp %in% group.ac)]

xx <- colorRampPalette(c("#fc8d59", "#ffffbf", "#91bfdb"))( 1+(max(group) - min(group)) )
group.cols <- xx[1+ (group - min(group))]

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
pdf(paste0("Figure_5_ML_longevity_COUNT.pdf"), width = 12, height = 6)
  layout(matrix(1:2, ncol = 2), width = c(2,0.3), height = c(1,1))
  par(mai=c(1, 1, 1, 0))
  hist(mRa, breaks = seq(min(mRa)-0.5, (max(mRa)+0.5), 1), col = group.cols, las = 2, xaxt = 'n',
       xlab = "# of (Increasing - Decreasing) alleles", 
       main = "Lifespan versus Haplotype Count")
  axis(1, at = seq(-10,10, 5), seq(-10,10, 5))

  par(mai=c(1, 0.2, 1, 0.2))
  legend_image <- as.raster(matrix(rev(xx), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
  text(x=1.5, y = seq(0, 1, l=5), labels = round(seq(min(group),max(group),l=5),0))
  rasterImage(legend_image, 0, 0, 1,1)
dev.off()

#plot(pull.pheno(mcross)[iix, "longevity"] ~ mRp2[iix])


pdf(paste0("Figure_5_ML_longevity_EFFECT.pdf"), width = 12, height = 6)
  plot(mRp2, observed, pch=19, cex=0.4, col = c("hotpink", "blue")[sex + 1], main = "Lifespan versus Haplotype Effect", xlab = "Estimated haplotype effect", ylab = "Lifespan", yaxt = "n")
  abline(lm(observed[sex==0] ~ mRp2[sex==0]), col = "hotpink", lwd=2, untf=T)
  abline(lm(observed[sex==1] ~ mRp2[sex==1]), col = "blue", lwd=2, untf=T)
  cor(observed[sex==0], mRp2[sex==0])
  cor(observed[sex==1], mRp2[sex==1])
  axis(2, seq(365, 1500, 4*15), seq(365, 1500, 4*15), las=2)
  legend("topleft", c("Female", "Male"), col = c("hotpink", "blue"), pch=19)
dev.off()


pdf(paste0("Figure_5_ML_longevity_EFFECT.pdf"), width = 28, height = 10)
par(cex=2)
par(cex.axis=1.5)
op <- par(cex=4)
layout(matrix(c(1,1,2,2,2,1,1,2,2,2), 2, 5, byrow = TRUE))

plot(c(1, 13), c(750, 1000), t = 'n', xaxt='n', xlab = "Allele Effect", ylab="Lifespan", las=2, yaxs='i')
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


width = 10

mRp2c <- cut(mRp2, breaks = seq(-150, 150, width), labels = FALSE)
iix <- which(pull.pheno(mcross)[, "longevity"] > 365)
aa <- boxplot(pull.pheno(mcross)[iix, "longevity"] ~ mRp2c[iix], plot=FALSE)
r <- 7:24#which(aa$n > 100)

xx <- colorRampPalette(c("#E41A1C", "white", "#377EB8"))( length(r) )

plot(c(1, length(r)), c(750, 1000), t = 'n', xaxt='n', xlab = "Expected effect on longevity ", ylab="Lifespan", las=2, yaxs='i')

phe <- pull.pheno(mcross)[iix, "longevity"]
sds <- c()
for(x in r){ sds <- c(sds, sd(phe[which(mRp2c[iix] == x)]) / sqrt(length(phe[which(mRp2c[iix] == x)]))); }


rect(1:length(r) - 0.25, rep(750,length(r)), 1:length(r) + 0.25, aa$stats[3,r], col = xx)
#points(aa$stats[3,], pch = 19, t = 'b')
#points(aa$stats[3,] + sds, pch = "-", cex=2, t = "h", lwd=2)
#points(aa$stats[3,] + sds, pch = "-", cex=3)
#points(aa$stats[3,] - sds, pch = "-", cex=2, t = "h", col = xx, lwd=2)
#points(aa$stats[3,] - sds, pch = "-", cex=3)
points(aa$conf[2,r], pch = "-", cex=2, t = "h", lwd=2)
points(aa$conf[2,r], pch = "-", cex=3)
points(aa$conf[1,r], pch = "-", cex=2, t = "h", col = xx, lwd=4)
points(aa$conf[1,r], pch = "-", cex=3)
mRp2c <- cut(mRp2, breaks = seq(-150, 150, width))

axis(1, at = 1:length(r), names(table(mRp2c))[min(r):max(r)])

dev.off()


