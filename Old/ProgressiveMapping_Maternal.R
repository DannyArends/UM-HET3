setwd("C:/Users/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Users/rqdt9/OneDrive - Northumbria University - Production Azure AD/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
gtsp <- pull.geno(mcross)

pORm <- names(which(lapply(lapply(apply(gtsp,2,unique,na.rm = TRUE), na.omit), length) == 2))
gtsp <- gtsp[, pORm]

matM <- names(which(apply(gtsp,2, function(x){any(5 %in% x)})))
gtsp <- gtsp[, matM]
gtsp[gtsp == 5] <- "A?"
gtsp[gtsp == 6] <- "B?"

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(gtsp), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(gtsp), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(gtsp)

# Our Progressive Mapping Sequence
msequence <- seq(365, 1100, 15)

lods.cM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  # Subset data based on the current threshold
  idx <- which(cdata[, "longevity"] >= x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  for(marker in colnames(gtsp)){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    idxF <- which(!is.na(mp))
    # Subset data based on the current marker
    sdata <- cdata[idxF,]
    mp <- mp[idxF]
    tryCatch(
      {
      lm.null <- lm(longevity ~ sex + site + cohort + treatment + 0, data = sdata)
      lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp + 0, data = sdata)
      n <- sum(!is.na(lm.alt$resid))
      lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
      lods.c <- c(lods.c, lod)
      },error=function(cond) {
        lods.c <<- c(lods.c, NA)
      }
    )
  }
  lods.cM <- rbind(lods.cM, lods.c)
}
colnames(lods.cM) <- colnames(gtsp)

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "Purples")

layout(matrix(c(1,1,1,1,1,1,2), ncol=7))
op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.cM), xlab="Maternal map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females and males", breaks = c(0, 2, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(">", minAge, " days")
axis(2, at = (1:length(yaxisD))[seq(1,length(yaxisD), 4)], yaxisD[seq(1,length(yaxisD), 4)], las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19))
box()
op <- par(mar = c(4,0,3,0))
plot(1:10, t = 'n', bty="n", xaxt ='n', yaxt ='n',xlab="")
legend("top", legend= c("<2.00", "2.00 - 3.65", "3.65 - 4.25", "4.25 - 4.95", "4.95 - 5.95", ">5.95"), fill = c("white", colz), bg="white")

### females

lods.fM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
  # Subset data based on the current threshold
  idx <- which(cdata[, "longevity"] >= x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  # Subset data to contain only females
  females <- which(cdata[, "sex"] == 0)
  cdata <- cdata[females,]
  gtsM <- gtsM[females, ]
  
  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  for(marker in colnames(gtsp)){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    idxF <- which(!is.na(mp))

    # Subset data based on the current marker
    sdata <- cdata[idxF,]
    mp <- mp[idxF]
    tryCatch(
      {
      lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = sdata)
      lm.alt <- lm(longevity ~ site + cohort + treatment + mp + 0, data = sdata)
      n <- sum(!is.na(lm.alt$resid))
      lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
      lods.c <- c(lods.c, lod)
      },error=function(cond) {
        lods.c <<- c(lods.c, NA)
      }
    )
  }
  lods.fM <- rbind(lods.fM, lods.c)
}
colnames(lods.fM) <- colnames(gtsp)

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "PuRd")[-1]

layout(matrix(c(1,1,1,1,1,1,2), ncol=7))
op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.fM), 1:nrow(lods.fM), t(lods.fM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females", breaks = c(0, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(">", minAge, " days")
axis(2, at = (1:length(yaxisD))[seq(1,length(yaxisD), 4)], yaxisD[seq(1,length(yaxisD), 4)], las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19))
box()
op <- par(mar = c(4,0,3,0))
plot(1:10, t = 'n', bty="n", xaxt ='n', yaxt ='n',xlab="")
legend("top", legend= c("<3.65", "3.65 - 4.25", "4.25 - 4.95", "4.95 - 5.95", ">5.95"), fill = c("white", colz), bg="white")

### males
lods.mM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] > x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  males <- which(cdata[, "sex"] == 1)
  cdata <- cdata[males,]
  gtsM <- gtsM[males, ]
  
  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  for(marker in colnames(gtsp)){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    idxF <- which(!is.na(mp))

    # Subset data based on the current marker
    sdata <- cdata[idxF,]
    mp <- mp[idxF]
    tryCatch(
      {
      lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = sdata)
      lm.alt <- lm(longevity ~ site + cohort + treatment + mp + 0, data = sdata)
      n <- sum(!is.na(lm.alt$resid))
      lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
      lods.c <- c(lods.c, lod)
      },error=function(cond) {
        lods.c <<- c(lods.c, NA)
      }
    )
  }
  lods.mM <- rbind(lods.mM, lods.c)
}
colnames(lods.mM) <- colnames(gtsp)

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "Blues")

layout(matrix(c(1,1,1,1,1,1,2), ncol=7))
op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.mM), 1:nrow(lods.mM), t(lods.mM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - males", breaks = c(0,2, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(">", minAge, " days")
axis(2, at = (1:length(yaxisD))[seq(1,length(yaxisD), 4)], yaxisD[seq(1,length(yaxisD), 4)], las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19))
box()
op <- par(mar = c(4,0,3,0))
plot(1:10, t = 'n', bty="n", xaxt ='n', yaxt ='n',xlab="")
legend("top", legend= c("<3.65", "3.65 - 4.25", "4.25 - 4.95", "4.95 - 5.95", ">5.95"), fill = c("white", colz), bg="white")

rownames(lods.mM) <- paste0("> ", msequence)
rownames(lods.fM) <- paste0("> ", msequence)
rownames(lods.cM) <- paste0("> ", msequence)

write.table(round(lods.mM,2), "progressiveMapping_mat_males.txt", sep = "\t", quote=FALSE)
write.table(round(lods.fM,2), "progressiveMapping_mat_females.txt", sep = "\t", quote=FALSE)
write.table(round(lods.cM,2), "progressiveMapping_mat_all.txt", sep = "\t", quote=FALSE)

threshold <- 3.65

mSign <- names(which(apply(lods.mM,2,function(x){any(x > threshold)})))
fSign <- names(which(apply(lods.fM,2,function(x){any(x > threshold)})))
cSign <- names(which(apply(lods.cM,2,function(x){any(x > threshold)})))

mSignificant <- colnames(lods.mM[which(apply(lods.mM[,mSign],2,max) > threshold), mSign])
fSignificant <- colnames(lods.fM[which(apply(lods.fM[,fSign],2,max) > threshold), fSign])
cSignificant <- colnames(lods.cM[which(apply(lods.cM[,cSign],2,max) > threshold), cSign])

#All
lapply(cSignificant, function(x){c(x, rownames(lods.cM)[which.max(lods.cM[,x])], max(lods.cM[,x]))})

lods.cM["> 650",]

#Females
lapply(fSignificant, function(x){c(x, rownames(lods.fM)[which.max(lods.fM[,x])], max(lods.fM[,x]))})

lods.fM["> 680",]

#Males
lapply(mSignificant, function(x){c(x, rownames(lods.mM)[which.max(lods.mM[,x])], max(lods.mM[,x]))})

lods.mM["> 605",]

getEffect <- function(mcross, gtsprob, timepoint = 365, sex = "all", marker = "1_3010274", model = "longevity ~ sex + site + cohort + treatment"){
  mp <- gtsprob[, grep(marker, colnames(gtsprob))]

  cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]),
                      geno = mp)
  idx <- which(cdata[, "longevity"] >= timepoint)
  cdata <- cdata[idx,]
  if(sex == 0 || sex == 1){ # males = 1, fem = 0
    cdata <- cdata[which(cdata[, "sex"] == sex),]
  }
  #cat("N=", nrow(cdata), "\n")

  gts <- cdata[,"geno"]
  mlm <- lm(as.formula(model), data = cdata)

  pheAdj <- rep(NA, nrow(cdata))
  names(pheAdj) <-  rownames(cdata)

  adj <- residuals(mlm) + mean(cdata[, "longevity"])
  pheAdj[names(adj)] <- adj
  OAmean <- mean(pheAdj[which(!is.na(gts))])
  means <- c(mean(pheAdj[which(gts == "A?")]),mean(pheAdj[which(gts == "B?")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "A?")]),std(pheAdj[which(gts == "B?")]))
  paste0(length(which(!is.na(gts))), "/",nrow(cdata), "\t", round(OAmean,0), "\t", paste0(round(means,0), " ± ", round(stderrs,2), collapse="\t"))
}

getEffect(mcross, gtsp, marker = "4_52524395", timepoint = 650)
getEffect(mcross, gtsp, marker = "4_74811205", timepoint = 590)
getEffect(mcross, gtsp, marker = "15_79013396", timepoint = 905)
getEffect(mcross, gtsp, marker = "17_18001459", timepoint = 605)
getEffect(mcross, gtsp, marker = "17_68770703", timepoint = 695)
getEffect(mcross, gtsp, marker = "18_52488251", timepoint = 365)

getEffect(mcross, gtsp, marker = "8_36994142", timepoint = 680, sex = 0, model = "longevity ~ site + cohort + treatment")
getEffect(mcross, gtsp, marker = "13_20905668", timepoint = 1100, sex = 0, model = "longevity ~ site + cohort + treatment")
getEffect(mcross, gtsp, marker = "18_52488251", timepoint = 650, sex = 0, model = "longevity ~ site + cohort + treatment")
getEffect(mcross, gtsp, marker = "18_67884842", timepoint = 380, sex = 0, model = "longevity ~ site + cohort + treatment")


getEffect(mcross, gtsp, marker = "4_52524395", timepoint = 605, sex = 1, model = "longevity ~ site + cohort + treatment")
getEffect(mcross, gtsp, marker = "6_54992703", timepoint = 1100, sex = 1, model = "longevity ~ site + cohort + treatment")
getEffect(mcross, gtsp, marker = "17_18001459", timepoint = 605, sex = 1, model = "longevity ~ site + cohort + treatment")
