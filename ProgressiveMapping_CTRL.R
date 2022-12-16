setwd("C:/Users/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Users/rqdt9/OneDrive - Northumbria University - Production Azure AD/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

# Our Progressive Mapping Sequence
msequence <- seq(365, 1100, 15)

lods.cM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      GenoID =pull.pheno(mcross)[, "GenoID"], 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] >= x & cdata[, "treatment"] == 0)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  lm.null <- lm(longevity ~ sex + site + cohort + 0, data = cdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(longevity ~ sex + site + cohort + mp + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.cM <- rbind(lods.cM, lods.c)
}
colnames(lods.cM) <- colnames(pull.geno(mcross))

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "Purples")

layout(matrix(c(1,1,1,1,1,1,2), ncol=7))
op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.cM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females and males", breaks = c(0, 2, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(">", minAge, " days")
axis(2, at = (1:length(yaxisD))[seq(1,length(yaxisD), 4)], yaxisD[seq(1,length(yaxisD), 4)], las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
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

  idx <- which(cdata[, "longevity"] >= x & cdata[, "treatment"] == 0)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  females <- which(cdata[, "sex"] == 0)
  cdata <- cdata[females,]
  gtsM <- gtsM[females, ]
  
  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  lm.null <- lm(longevity ~ site + cohort + 0, data = cdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(longevity ~ site + cohort + mp + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.fM <- rbind(lods.fM, lods.c)
}
colnames(lods.fM) <- colnames(pull.geno(mcross))

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
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
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

  idx <- which(cdata[, "longevity"] > x & cdata[, "treatment"] == 0)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  males <- which(cdata[, "sex"] == 1)
  cdata <- cdata[males,]
  gtsM <- gtsM[males, ]
  
  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  lm.null <- lm(longevity ~ site + cohort + 0, data = cdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(longevity ~ site + cohort + mp + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.mM <- rbind(lods.mM, lods.c)
}
colnames(lods.mM) <- colnames(pull.geno(mcross))

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
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
box()
op <- par(mar = c(4,0,3,0))
plot(1:10, t = 'n', bty="n", xaxt ='n', yaxt ='n',xlab="")
legend("top", legend= c("<3.65", "3.65 - 4.25", "4.25 - 4.95", "4.95 - 5.95", ">5.95"), fill = c("white", colz), bg="white")

rownames(lods.mM) <- paste0("> ", msequence)
rownames(lods.fM) <- paste0("> ", msequence)
rownames(lods.cM) <- paste0("> ", msequence)

write.table(round(lods.mM,2), "progressiveMapping_ctrl_males.txt", sep = "\t", quote=FALSE)
write.table(round(lods.fM,2), "progressiveMapping_ctrl_females.txt", sep = "\t", quote=FALSE)
write.table(round(lods.cM,2), "progressiveMapping_ctrl_all.txt", sep = "\t", quote=FALSE)

threshold <- 3.65

mSign <- names(which(apply(lods.mM,2,function(x){any(x > threshold)})))
fSign <- names(which(apply(lods.fM,2,function(x){any(x > threshold)})))
cSign <- names(which(apply(lods.cM,2,function(x){any(x > threshold)})))

mSignificant <- colnames(lods.mM[which(apply(lods.mM[,mSign],1,max) > threshold), mSign])
fSignificant <- colnames(lods.fM[which(apply(lods.fM[,fSign],1,max) > threshold), fSign])
cSignificant <- colnames(lods.cM[which(apply(lods.cM[,cSign],1,max) > threshold), cSign])

#All
lapply(cSignificant, function(x){c(x, rownames(lods.cM)[which.max(lods.cM[,x])], max(lods.cM[,x]))})

lods.cM["> 740",]

#Females
lapply(fSignificant, function(x){c(x, rownames(lods.fM)[which.max(lods.fM[,x])], max(lods.fM[,x]))})

round(lods.fM["> 500",],2)

#Males
lapply(mSignificant, function(x){c(x, rownames(lods.mM)[which.max(lods.mM[,x])], max(lods.mM[,x]))})

round(lods.mM["> 770",],2)

getEffect <- function(mcross, gtsprob, timepoint = 365, sex = "all", marker = "1_3010274", model = "longevity ~ sex + site + cohort"){
  mp <- gtsprob[, grep(marker, colnames(gtsprob))]
  gts <- unlist(lapply(lapply(lapply(apply(mp,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]),
                      geno = gts)
  idx <- which(cdata[, "longevity"] >= timepoint & cdata[, "treatment"] == 0)
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
  means <- c(mean(pheAdj[which(gts == "AC")]),mean(pheAdj[which(gts == "AD")]),mean(pheAdj[which(gts == "BC")]),mean(pheAdj[which(gts == "BD")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "AC")]),std(pheAdj[which(gts == "AD")]),std(pheAdj[which(gts == "BC")]),std(pheAdj[which(gts == "BD")]))
  paste0(length(which(!is.na(gts))), "/",nrow(cdata), "\t", round(OAmean,0), "\t", paste0(round(means,0), " ± ", round(stderrs,2), collapse="\t"))
}

getEffect(mcross, gtsp, marker = "1_58265070", timepoint = 740)
getEffect(mcross, gtsp, marker = "2_83915237", timepoint = 815)
getEffect(mcross, gtsp, marker = "3_103803615", timepoint = 770)
getEffect(mcross, gtsp, marker = "9_38672644", timepoint = 890)
getEffect(mcross, gtsp, marker = "9_83925296", timepoint = 1040)
getEffect(mcross, gtsp, marker = "12_117441622", timepoint = 530)
getEffect(mcross, gtsp, marker = "13_118386771", timepoint = 410)
getEffect(mcross, gtsp, marker = "14_79402836", timepoint = 605)
getEffect(mcross, gtsp, marker = "17_18001459", timepoint = 590)
getEffect(mcross, gtsp, marker = "19_17436579", timepoint = 1025)

getEffect(mcross, gtsp, marker = "2_148442635", timepoint = 500, sex = 0, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "3_103803615", timepoint = 770, sex = 0, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "9_29939029", timepoint = 470, sex = 0, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "16_50590634", timepoint = 365, sex = 0, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "19_17436579", timepoint = 1025, sex = 0, model = "longevity ~ site + cohort")

getEffect(mcross, gtsp, marker = "1_180168619", timepoint = 770, sex = 1, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "2_83915237", timepoint = 905, sex = 1, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "3_83838529", timepoint = 1070, sex = 1, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "4_149128877", timepoint = 650, sex = 1, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "9_124056586", timepoint = 815, sex = 1, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "10_72780332", timepoint = 995, sex = 1, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "12_112855820", timepoint = 860, sex = 1, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "14_20686566", timepoint = 1040, sex = 1, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "15_78849558", timepoint = 365, sex = 1, model = "longevity ~ site + cohort")
getEffect(mcross, gtsp, marker = "17_18001459", timepoint = 590, sex = 1, model = "longevity ~ site + cohort")




