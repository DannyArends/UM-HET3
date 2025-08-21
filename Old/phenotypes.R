#
# phenotypes.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Old pre-processing of phenotype data, with QTL scans to see if adjustment improves
#

setwd("C:/Github/UM-HET3/files/merged")

gts4way <- read.table("gts4way.Juli2021.txt", sep="\t",na.strings=c("","NA", "??", "XX"))
gts4wayRqtl <- read.table("gts4way.rqtl.Juli2021.txt", sep="\t")
map <- read.table("map.gts4way.Juli2021.txt", sep="\t")
ind <- read.table("ind.gts4way.Juli2021.txt", sep="\t")

map[map[, "Chr"] == "X", "Chr"] <- 20

gts4wayNum <- apply(gts4wayRqtl, 2, as.numeric)
rownames(gts4wayNum) <- rownames(gts4wayNum)

# Phenotype QC, remove individuals with missing start weight and/or no longevity data
hasBW <- rownames(ind)[which(!is.na(ind[, "BodyWeight_HET3_ITP_6m"]) & !is.na(ind[, "Longevity_HET3_ITP"]))]
ind <- ind[hasBW,]
gts4way <- gts4way[,hasBW]

dim(ind)
dim(gts4way)

# Phenotype QC, remove individuals with missing data or data out of range (mean +/- 2SD)
res <- c()
validInd <- c()
for(site in unique(ind[,"Site"])){
  malesOnSite <- ind[which(ind$Site == site & ind$Sex == "M"),]
  femalesOnSite <- ind[which(ind$Site == site & ind$Sex == "F"),]
  resRow <- c()
  for(pheCol in c(8, 9, 10, 11)){
    phe <- colnames(ind)[pheCol]
    cat(site, phe, "\n")
    mean.m <- round(mean(malesOnSite[,phe], na.rm = TRUE),1)
    sds.m <- round(sd(malesOnSite[,phe], na.rm = TRUE),1)
    isValid.m <- rownames(malesOnSite)[which((malesOnSite[,phe] > mean.m - (sds.m * 2)) & (malesOnSite[,phe] < mean.m + (sds.m * 2)))]
    cat("males:", length(isValid.m), "/", nrow(malesOnSite), " mean: ",mean.m, "+/-",sds.m, "\n")

    mean.f <- round(mean(femalesOnSite[,phe], na.rm = TRUE),1)
    sds.f <- round(sd(femalesOnSite[,phe], na.rm = TRUE),1)
    isValid.f <- rownames(femalesOnSite)[which((femalesOnSite[,phe] > mean.f - (sds.f * 2)) & (femalesOnSite[,phe] < mean.f + (sds.f * 2)))]
    cat("females:", length(isValid.f), "/", nrow(femalesOnSite), " mean: ",mean.f, "+/-",sds.f, "\n")

    resRow <- c(resRow, nrow(malesOnSite), mean.m, sds.m, nrow(femalesOnSite), mean.f, sds.f)
    validInd <- rbind(validInd, c(site, phe, paste0(isValid.m, collapse=","),paste0(isValid.f, collapse=",")))
  }
  res <- rbind(res, resRow)
}
rownames(res) <- unique(ind[,"Site"])
colnames(res) <- gsub("BodyWeight_HET3_ITP_", "BW", unlist(lapply(colnames(ind)[c(8:11)], paste, c("N.M", "Mean.M", "SD.M", "N.F", "Mean.F", "SD.F"))))

colnames(validInd) <- c("Site", "Phenotype", "Males", "Females")

validBW <- validInd[grep("BodyWeight", validInd[,2]),]
indAdjusted <- matrix(NA, nrow(ind), ncol(ind), dimnames = list(rownames(ind), colnames(ind)))
indAdjusted <- data.frame(indAdjusted)
indAdjusted[,1] <- ind[,1]
indAdjusted[,2] <- ind[,2]
indAdjusted[,3] <- ind[,3]
indAdjusted[,4] <- ind[,4]
indAdjusted[,6] <- ind[,6]
indAdjusted[,7] <- ind[,7]
indAdjusted <- indAdjusted[,-5]

for(pheCol in c(7:10)){
  phe <- colnames(indAdjusted)[pheCol]
  cat("Adjusting", phe, "\n")
  # Individuals valid for mapping
  
  isValid.M <- unlist(lapply(validInd[which(validInd[, "Phenotype"] == phe),"Males"], strsplit, ","))
  isValid.F <- unlist(lapply(validInd[which(validInd[, "Phenotype"] == phe),"Females"], strsplit, ","))
  mformula <- as.formula(paste0(phe, " ~ as.factor(Site) + as.factor(Treatment_Effect) + as.factor(Cohort.Year)"))
  
  indAdjusted[isValid.M, phe] <- as.numeric(ind[isValid.M, phe])
  adj.M <- lm(mformula, indAdjusted[isValid.M,])  # Adjust phenotypes for males
  BWAdj.M <- coef(adj.M)["(Intercept)"] + residuals(adj.M)
  indAdjusted[names(BWAdj.M), phe] <- round(BWAdj.M,2)
  
  indAdjusted[isValid.F, phe] <- as.numeric(ind[isValid.F, phe])
  adj.F <- lm(mformula, indAdjusted[isValid.F,])  # Adjust phenotypes for females
  BWAdj.F  <- coef(adj.F)["(Intercept)"] + residuals(adj.F)
  indAdjusted[names(BWAdj.F), phe] <- round(BWAdj.F,2)
}

write.table(indAdjusted, file = "bodyweights.nooutliers.adjusted.txt", quote=FALSE, sep="\t", na = "")


# Initial QTL mapping (I adjust again for some reason)
lods <- c()
for(pheCol in c(8, 9, 10, 11)){
  phe <- colnames(ind)[pheCol]

  # Individuals valid for mapping
  isValid.M <- unlist(lapply(validInd[which(validInd[, "Phenotype"] == phe),"Males"], strsplit, ","))
  isValid.F <- unlist(lapply(validInd[which(validInd[, "Phenotype"] == phe),"Females"], strsplit, ","))

  # Model (Site, Treatment, Cohort)
  mformula <- as.formula(paste0(phe, " ~ as.factor(Site) + as.factor(Treatment_Effect) + as.factor(Cohort.Year)"))

  adj.M <- lm(mformula, ind[isValid.M,])  # Adjust phenotypes for males
  BWAdj.M <- coef(adj.M)["(Intercept)"] + residuals(adj.M)

  adj.F <- lm(mformula, ind[isValid.F,])  # Adjust phenotypes for females
  BWAdj.F  <- coef(adj.F)["(Intercept)"] + residuals(adj.F)

  # Map and plot males all markers
  op <- par(mfrow=c(2, 1))
  mm <- apply(gts4way[,names(BWAdj.M)], 1, function(gts){ 
    if(length(table(gts)) > 1) return(anova(lm(BWAdj.M ~ gts))[[5]][1])
    return(NA)
  })
  plot(-log10(mm), t = 'b', main = paste0("Male ",phe,"(Adj: Site & Cohort)"))
  lods <- rbind(lods, -log10(mm))

  # Map and plot females all markers
  mm <- apply(gts4way[,names(BWAdj.F)], 1, function(gts){ 
    if(length(table(gts)) > 1) return(anova(lm(BWAdj.F ~ gts))[[5]][1])
    return(NA)
  })
  plot(-log10(mm), t = 'b', main = paste0("Female ",phe,"(Adj: Site & Cohort)"))
  lods <- rbind(lods, -log10(mm))
}
rownames(lods) <- gsub("BodyWeight_HET3_ITP_", "BW", unlist(lapply(colnames(ind)[c(8:11)], paste, c("M", "F"))))

# Paternal and Maternal maps
mMap <- rownames(map)[grep("Mat", map[,"type"])]
fMap <- rownames(map)[grep("Pat", map[,"type"])]

op <- par(mfrow=c(2,1))

library(RColorBrewer)
mcol = brewer.pal(4, "PuRd")

# Plot paternal map
op <- par(mar = c(4,10,3,1))
op <- par(cex = 1)
image(1:length(fMap), 1:8, t(lods[c(1, 3, 5, 7, 2, 4, 6, 8),fMap]), xlab="Paternal Map",yaxt='n',ylab="",xaxt='n')
abline(v=which(diff(as.numeric(map[fMap,"Chr"])) != 0))
abline(h=c(4.5), lty=2)
axis(2, at = 1:8, rownames(lods)[c(1, 3, 5, 7, 2, 4, 6, 8)], las = 2)
chrP <- c(0, which(diff(as.numeric(map[fMap,"Chr"])) != 0))
axis(1, at = chrP + (diff(c(chrP, length(fMap)))/2), 1:19)
box()

# Plot maternal map
op <- par(mar = c(4,10,3,1))
image(1:length(mMap), 1:8, t(lods[c(1, 3, 5, 7, 2, 4, 6, 8),mMap]), xlab="Maternal Map",yaxt='n',ylab="",xaxt='n')
abline(v=which(diff(as.numeric(map[mMap,"Chr"])) != 0))
abline(h=c(4.5), lty=2)
axis(2, at = 1:8, rownames(lods)[c(1, 3, 5, 7, 2, 4, 6, 8)], las = 2)
chrP <- c(0, which(diff(as.numeric(map[mMap,"Chr"])) != 0))
axis(1, at = chrP + (diff(c(chrP, length(mMap)))/2), 1:19)
box()

### QTL mapping of the longevity phenotypes in Males
phe <- "Longevity_HET3_ITP"
males <- ind[which(ind[,"Sex"] == "M"),]

adj.M <- lm(Longevity_HET3_ITP ~ as.factor(Site) + as.factor(Treatment_Effect) + as.factor(Cohort.Year), males)
phe.M <- residuals(adj.M) + coef(adj.M)["(Intercept)"]

lods.M <- c()
for(todrop in seq(0, 0.6, 0.1)){
  sorted <- sort(phe.M)
  phe.M.dropped <- sorted[round(length(sorted) * todrop):length(sorted)]
  cat("Left with:", length(phe.M.dropped), "\n")
  # Map the remaining males
  mm <- apply(gts4way[,names(phe.M.dropped)], 1, function(gts){ 
    if(length(table(gts)) > 1) return(anova(lm(phe.M.dropped ~ gts))[[5]][1])
    return(NA)
  })
  lods.M <- rbind(lods.M, -log10(mm))
}

op <- par(mfrow=c(2,1))

# Plot maternal map
op <- par(mar = c(4,6,3,1))
image(1:length(mMap), 1:7, t(lods.M[,mMap]), xlab="",yaxt='n',ylab="",xaxt='n')
abline(v=which(diff(as.numeric(map[mMap,"Chr"])) != 0))
axis(2, at = 1:7, paste0("Top ", 100 - (seq(0, 0.6, 0.1) * 100), "%"), las = 2)
chrP <- c(0, which(diff(as.numeric(map[mMap,"Chr"])) != 0))
axis(1, at = chrP + (diff(c(chrP, length(mMap)))/2), 1:19)
box()

# Plot paternal map
op <- par(mar = c(4,6,3,1))
op <- par(cex = 1)
image(1:length(fMap), 1:7, t(lods.M[,fMap]), xlab="",yaxt='n',ylab="",xaxt='n')
abline(v=which(diff(as.numeric(map[fMap,"Chr"])) != 0))
axis(2, at = 1:7, paste0("Top ", 100 - (seq(0, 0.6, 0.1) * 100), "%"), las = 2)
chrP <- c(0, which(diff(as.numeric(map[fMap,"Chr"])) != 0))
axis(1, at = chrP + (diff(c(chrP, length(fMap)))/2), 1:19)
box()


### QTL mapping of the longevity phenotypes in Males
females <- ind[which(ind[,"Sex"] == "F"),]

adj.F <- lm(Longevity_HET3_ITP ~ as.factor(Site) + as.factor(Treatment_Effect) + as.factor(Cohort.Year), females)
phe.F <- residuals(adj.F) + coef(adj.F)["(Intercept)"]

lods.F <- c()
for(todrop in seq(0, 0.6, 0.1)){
  sorted <- sort(phe.F)
  phe.F.dropped <- sorted[round(length(sorted) * todrop):length(sorted)]
  cat("Left with:", length(phe.F.dropped), "\n")
  # Map the remaining males
  mm <- apply(gts4way[,names(phe.F.dropped)], 1, function(gts){ 
    if(length(table(gts)) > 1) return(anova(lm(phe.F.dropped ~ gts))[[5]][1])
    return(NA)
  })
  lods.F <- rbind(lods.F, -log10(mm))
}

op <- par(mfrow=c(2,1))

# Plot maternal map
op <- par(mar = c(4,6,3,1))
image(1:length(mMap), 1:7, t(lods.F[,mMap]), xlab="",yaxt='n',ylab="",xaxt='n')
abline(v=which(diff(as.numeric(map[mMap,"Chr"])) != 0))
axis(2, at = 1:7, paste0("Top ", 100 - (seq(0, 0.6, 0.1) * 100), "%"), las = 2)
chrP <- c(0, which(diff(as.numeric(map[mMap,"Chr"])) != 0))
axis(1, at = chrP + (diff(c(chrP, length(mMap)))/2), 1:19)
box()

# Plot paternal map
op <- par(mar = c(4,6,3,1))
op <- par(cex = 1)
image(1:length(fMap), 1:7, t(lods.F[,fMap]), xlab="",yaxt='n',ylab="",xaxt='n')
abline(v=which(diff(as.numeric(map[fMap,"Chr"])) != 0))
axis(2, at = 1:7, paste0("Top ", 100 - (seq(0, 0.6, 0.1) * 100), "%"), las = 2)
chrP <- c(0, which(diff(as.numeric(map[fMap,"Chr"])) != 0))
axis(1, at = chrP + (diff(c(chrP, length(fMap)))/2), 1:19)
box()

