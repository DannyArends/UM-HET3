setwd("C:/Github/UM-HET3/files/merged")

gts4way <- read.table("gts4way.txt", sep="\t",na.strings=c("","NA", "??", "XX"))
gts4wayRqtl <- read.table("gts4way.rqtl.txt", sep="\t")
map <- read.table("map.gts4way.txt", sep="\t")
ind <- read.table("ind.gts4way.txt", sep="\t")

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
  for(pheCol in c(5, 8, 9, 10, 11)){
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
colnames(res) <- gsub("BodyWeight_HET3_ITP_", "BW", unlist(lapply(colnames(ind)[c(5,8:11)], paste, c("N.M", "Mean.M", "SD.M", "N.F", "Mean.F", "SD.F"))))

colnames(validInd) <- c("Site", "Phenotype", "Males", "Females")

# Initial QTL mapping
lods <- c()
for(pheCol in c(5, 8, 9, 10, 11)){
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
  plot(-log10(mm), col=map[names(mm), "Chr"], t = 'b', main = paste0("Male ",phe,"(Adj: Site & Cohort)"))
  lods <- rbind(lods, -log10(mm))

  # Map and plot females all markers
  mm <- apply(gts4way[,names(BWAdj.F)], 1, function(gts){ 
    if(length(table(gts)) > 1) return(anova(lm(BWAdj.F ~ gts))[[5]][1])
    return(NA)
  })
  plot(-log10(mm), col=map[names(mm), "Chr"], t = 'b', main = paste0("Female ",phe,"(Adj: Site & Cohort)"))
  lods <- rbind(lods, -log10(mm))
}
rownames(lods) <- gsub("BodyWeight_HET3_ITP_", "BW", unlist(lapply(colnames(ind)[c(5,8:11)], paste, c("M", "F"))))

# Paternal and Maternal maps
mMap <- rownames(map)[grep("Mat", map[,"type"])]
fMap <- rownames(map)[grep("Pat", map[,"type"])]

op <- par(mfrow=c(2,1))
# Plot paternal map
op <- par(mar = c(4,10,3,1))
op <- par(cex = 0.5)
image(1:length(fMap), 1:10, t(lods[c(1, 2, 3, 5, 7, 9, 4, 6, 8, 10),fMap]), xlab="Paternal Map",yaxt='n',ylab="",xaxt='n')
abline(v=which(diff(map[fMap,"Chr"]) != 0))
abline(h=c(2.5, 6.5), lty=2)
axis(2, at = 1:10, rownames(lods)[c(1, 2, 3, 5, 7, 9, 4, 6, 8, 10)], las = 2)
chrP <- c(0, which(diff(map[fMap,"Chr"]) != 0))
axis(1, at = chrP + (diff(c(chrP, length(fMap)))/2), 1:19)
box()

# Plot maternal map
op <- par(mar = c(4,10,3,1))
image(1:length(mMap), 1:10, t(lods[c(1, 2, 3, 5, 7, 9, 4, 6, 8, 10),mMap]), xlab="Maternal Map",yaxt='n',ylab="",xaxt='n')
abline(v=which(diff(map[mMap,"Chr"]) != 0))
abline(h=c(2.5, 6.5), lty=2)
axis(2, at = 1:10, rownames(lods)[c(1, 2, 3, 5, 7, 9, 4, 6, 8, 10)], las = 2)
chrP <- c(0, which(diff(map[mMap,"Chr"]) != 0))
axis(1, at = chrP + (diff(c(chrP, length(mMap)))/2), 1:19)
box()

plot(x = c(0, 20), y = c(0, max(map[, "Position"])), t = 'n', xlab="Chromosome", ylab="Pos")
lapply(rev(colnames(lods)), function(mname){
  isMat = (2*grepl("Mat", map[mname, "type"])) - 1
  mType = c("Green", "Blue", "Gold", "red")[as.numeric(factor(map[mname, "type"], levels = c("Mat1", "Mat2", "Pat1", "Pat2")))]
  rect(map[mname, "Chr"] + isMat * 0, 0, map[mname, "Chr"] + isMat * 0.2, map[mname,"Position"], col= mType, border="white")
})


lapply(colnames(lods), function(mname){
  isMat = (2*grepl("Mat", map[mname, "type"])) - 1
  points(c(map[mname, "Chr"] + isMat * 0.1,map[mname, "Chr"] + isMat * 0.1), c(map[mname, "Position"], map[mname, "Position"]), pch = '-', cex= as.numeric(lods[2, mname]))
})
