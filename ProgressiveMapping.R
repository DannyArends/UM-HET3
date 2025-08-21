#
# ProgressiveMapping.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Map longevity QTLs by progressively raising the threshold for inclusion (T-age) by 15 days (starting at 20 days)
# We perform a full scan at each T-age for both Combined (F+M)), Females, and Males
# We write out the results in a matrix , and visualize
#

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

# Writing out 
iix <- which(pull.pheno(mcross)[,"longevity"] >= 365)
write.table(pull.pheno(mcross)[iix,"GenoID"], "Cases_UM_HET3.txt",sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

write.table(map, "genetic_map.txt", sep = "\t", quote=FALSE)

# Our Progressive Mapping Sequence
msequence <- seq(20, 1100, 15)
markers <- unique(unlist(lapply(strsplit(colnames(gtsp), ":"), "[",1)))

lods.cM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.factor(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] >= x & cdata[, "sex"] == 1)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  lm.null <- lm(longevity ~ sex + site + cohort + treatment + 0, data = cdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.cM <- rbind(lods.cM, lods.c)
  cat("Done", x, "\n")
}
colnames(lods.cM) <- colnames(pull.geno(mcross))
rownames(lods.cM) <- paste0("> ", msequence)

write.table(round(lods.cM,2), "progressiveMapping_all.txt", sep = "\t", quote=FALSE)

subset <- map[which(as.numeric(map[,1]) %in% 1:4),]
subset <- cbind(subset, cpos = NA)
gap <- 10000000
chr.start <- c(0)
chr.mids <- c()
cp <- 0
for(x in 1:4){
  cl <- max(as.numeric(subset[which(subset[,1] == x),2]))
  chr.start <- c(chr.start, cl + cp + gap)
  subset[which(subset[,1] == x), "cpos"] <- as.numeric(subset[which(subset[,1] == x), 2]) + cp
  chr.mids <- c(chr.mids, chr.start[x] + cl/2)
  cp = cl + cp + gap
}

plot(c(0, max(chr.start)), y = c(0, 7), t = 'n', ylab = "LOD", xlab = "Chr",xaxt="n", las=2)
for(x in 1:4){
  points(subset[which(subset[,1] == x),"cpos"], lods.cM[, rownames(subset)[which(subset[,1] == x)]], t = 'l')
}

axis(1, at = chr.mids, paste0("Chr", 1:4))

png("ProgressiveMapping_UMHET3.png", width = 1400, height = 1050)

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "Purples")

layout(matrix(c(1,1,1,1,1,1,2,3,3,3,3,3,3,4,5,5,5,5,5,5,6), ncol=7, byrow=TRUE))
op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.cM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females and males", breaks = c(0, 2, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
#contour(1:ncol(lods.cM), 1:nrow(lods.cM), t(lods.cM), levels = c(0, 3.65, 4.25, 4.95, 5.95, 100), add = TRUE)

abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(">", msequence, " days")
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

  idx <- which(cdata[, "longevity"] >= x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  females <- which(cdata[, "sex"] == 0)
  cdata <- cdata[females,]
  gtsM <- gtsM[females, ]
  
  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = cdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(longevity ~ site + cohort + treatment + mp + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.fM <- rbind(lods.fM, lods.c)
  cat("Done", x, "\n")
}
colnames(lods.fM) <- colnames(pull.geno(mcross))

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "PuRd")[-1]

op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.fM), 1:nrow(lods.fM), t(lods.fM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - females", breaks = c(0, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(">", msequence, " days")
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

  idx <- which(cdata[, "longevity"] >= x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  males <- which(cdata[, "sex"] == 1)
  cdata <- cdata[males,]
  gtsM <- gtsM[males, ]
  
  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = cdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(longevity ~ site + cohort + treatment + mp + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.mM <- rbind(lods.mM, lods.c)
  cat("Done", x, "\n")
}
colnames(lods.mM) <- colnames(pull.geno(mcross))

# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "Blues")

op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.mM), 1:nrow(lods.mM), t(lods.mM), xlab="4-way map",
      yaxt='n',ylab="",xaxt='n', main="Longevity - males", breaks = c(0,2, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))
abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(">", msequence, " days")
axis(2, at = (1:length(yaxisD))[seq(1,length(yaxisD), 4)], yaxisD[seq(1,length(yaxisD), 4)], las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
box()
op <- par(mar = c(4,0,3,0))
plot(1:10, t = 'n', bty="n", xaxt ='n', yaxt ='n',xlab="")
legend("top", legend= c("<3.65", "3.65 - 4.25", "4.25 - 4.95", "4.95 - 5.95", ">5.95"), fill = c("white", colz), bg="white")

dev.off()


rownames(lods.mM) <- paste0("> ", msequence)
rownames(lods.fM) <- paste0("> ", msequence)
rownames(lods.cM) <- paste0("> ", msequence)

write.table(round(lods.mM,2), "progressiveMapping_males20D.txt", sep = "\t", quote=FALSE)
write.table(round(lods.fM,2), "progressiveMapping_females20D.txt", sep = "\t", quote=FALSE)
write.table(round(lods.cM,2), "progressiveMapping_all20D.txt", sep = "\t", quote=FALSE)

threshold <- 3.65

mSign <- names(which(apply(lods.mM,2,function(x){any(x > threshold)})))
fSign <- names(which(apply(lods.fM,2,function(x){any(x > threshold)})))
cSign <- names(which(apply(lods.cM,2,function(x){any(x > threshold)})))

mSignificant <- colnames(lods.mM[which(apply(lods.mM[,mSign],2,max) > threshold), mSign])
fSignificant <- colnames(lods.fM[which(apply(lods.fM[,fSign],2,max) > threshold), fSign])
cSignificant <- colnames(lods.cM[which(apply(lods.cM[,cSign],2,max) > threshold), cSign])

#All
lapply(cSignificant, function(x){c(x, rownames(lods.cM)[which.max(lods.cM[,x])], max(lods.cM[,x]))})

lods.cM["> 860",]

#Females
lapply(fSignificant, function(x){c(x, rownames(lods.fM)[which.max(lods.fM[,x])], max(lods.fM[,x]))})

lods.fM["> 770",]

#Males
lapply(mSignificant, function(x){c(x, rownames(lods.mM)[which.max(lods.mM[,x])], max(lods.mM[,x]))})

lods.mM["> 770",]

getEffect <- function(mcross, gtsprob, timepoint = 365, sex = "all", marker = "1_3010274", model = "longevity ~ sex + site + cohort + treatment"){
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
  means <- c(mean(pheAdj[which(gts == "AC")]),mean(pheAdj[which(gts == "AD")]),mean(pheAdj[which(gts == "BC")]),mean(pheAdj[which(gts == "BD")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "AC")]),std(pheAdj[which(gts == "AD")]),std(pheAdj[which(gts == "BC")]),std(pheAdj[which(gts == "BD")]))
  paste0(length(which(!is.na(gts))), "/",nrow(cdata), "\t", round(OAmean,1), "\t", paste0(round(means,1), " ± ", round(stderrs,2), collapse="\t"))
}

getEffect(mcross, gtsp, marker = "1_3010272", timepoint = 860)        # Vita1A
getEffect(mcross, gtsp, marker = "1_24042124", timepoint = 695)       # Vita1B
getEffect(mcross, gtsp, marker = "9_104091597", timepoint = 1025)     # Vita9b
getEffect(mcross, gtsp, marker = "9_124056586", timepoint = 980)      # Vita9c
getEffect(mcross, gtsp, marker = "10_72780332", timepoint = 980)      # Vita10a
getEffect(mcross, gtsp, marker = "12_112855820", timepoint = 635)     # Vita12a
getEffect(mcross, gtsp, marker = "13_83858506", timepoint = 230)      # Vita13a
getEffect(mcross, gtsp, marker = "14_78415875", timepoint = 35)       # Vita14a
getEffect(mcross, gtsp, marker = "14_101437466", timepoint = 860)     # Vita14b
getEffect(mcross, gtsp, marker = "15_74248242", timepoint = 905)      # Vita15a
getEffect(mcross, gtsp, marker = "17_32883804", timepoint = 665)      # Vita17a
getEffect(mcross, gtsp, marker = "18_52488251", timepoint = 140)      # Vita18a

#Female
getEffect(mcross, gtsp, marker = "2_112255823", timepoint = 800, sex = 0, model = "longevity ~ site + cohort + treatment")  # Vita2b
getEffect(mcross, gtsp, marker = "2_148442635", timepoint = 545, sex = 0, model = "longevity ~ site + cohort + treatment")  # Vita2c
getEffect(mcross, gtsp, marker = "9_34932404", timepoint = 95, sex = 0, model = "longevity ~ site + cohort + treatment")    # Vita9a
getEffect(mcross, gtsp, marker = "11_82176894", timepoint = 1040, sex = 0, model = "longevity ~ site + cohort + treatment") # Vita11b
getEffect(mcross, gtsp, marker = "X_156343080", timepoint = 740, sex = 0, model = "longevity ~ site + cohort + treatment")  # VitaXb

#Males
getEffect(mcross, gtsp, marker = "1_121483290", timepoint = 230, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita1c
getEffect(mcross, gtsp, marker = "1_167148678", timepoint = 230, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita1d
getEffect(mcross, gtsp, marker = "2_89156987", timepoint = 140, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita2a
getEffect(mcross, gtsp, marker = "3_83354281", timepoint = 170, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita3a
getEffect(mcross, gtsp, marker = "4_52524395", timepoint = 650, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita4a
getEffect(mcross, gtsp, marker = "4_154254581", timepoint = 630, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita4b
getEffect(mcross, gtsp, marker = "5_67573068", timepoint = 1085, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita5a
getEffect(mcross, gtsp, marker = "6_93680853", timepoint = 140, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita6a
getEffect(mcross, gtsp, marker = "6_132762500", timepoint = 320, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita6b
getEffect(mcross, gtsp, marker = "11_6599922", timepoint = 635, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita11a
getEffect(mcross, gtsp, marker = "11_113729074", timepoint = 35, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita11c
getEffect(mcross, gtsp, marker = "15_99306167", timepoint = 275, sex = 1, model = "longevity ~ site + cohort + treatment") # Vita15b
getEffect(mcross, gtsp, marker = "X_36008085", timepoint = 260, sex = 1, model = "longevity ~ site + cohort + treatment") # VitaXa











