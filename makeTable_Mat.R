library(RColorBrewer)
library(svglite)

setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
map <- read.table("genetic_map.txt", sep = "\t")

lods.mM <- read.table("progressiveMapping_mat_males.txt", sep = "\t", check.names=FALSE)
lods.fM <- read.table("progressiveMapping_mat_females.txt", sep = "\t", check.names=FALSE)
lods.cM <- read.table("progressiveMapping_mat_all.txt", sep = "\t", check.names=FALSE)

threshold <- 3.65

mSign <- names(which(apply(lods.mM,2,function(x){any(x > threshold)})))
fSign <- names(which(apply(lods.fM,2,function(x){any(x > threshold)})))
cSign <- names(which(apply(lods.cM,2,function(x){any(x > threshold)})))

mSignificant <- colnames(lods.mM[which(apply(lods.mM[,mSign],2,max) > threshold), mSign])
fSignificant <- colnames(lods.fM[which(apply(lods.fM[,fSign],2,max) > threshold), fSign])
cSignificant <- colnames(lods.cM[which(apply(lods.cM[,cSign],2,max) > threshold), cSign])


lapply(cSignificant, function(x){c(x, rownames(lods.cM)[which.max(lods.cM[,x])], max(lods.cM[,x]))})
lapply(fSignificant, function(x){c(x, rownames(lods.fM)[which.max(lods.fM[,x])], max(lods.fM[,x]))})
lapply(mSignificant, function(x){c(x, rownames(lods.mM)[which.max(lods.mM[,x])], max(lods.mM[,x]))})


setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

getEffect <- function(mcross, gtsprob, timepoint = 365, sex = "all", marker = "1_3010274", model = "longevity ~ sex + site + cohort + treatment"){
  mp <- gtsprob[, grep(marker, colnames(gtsprob))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  gts[gts == "AC"] <- "A"; gts[gts == "AD"] <- "A"
  gts[gts == "BC"] <- "B"; gts[gts == "BD"] <- "B"

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
  means <- c(mean(pheAdj[which(gts == "A")]),mean(pheAdj[which(gts == "B")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "A")]),std(pheAdj[which(gts == "B")]))
  paste0(length(which(!is.na(gts))), "/",nrow(cdata), "\t", round(OAmean,0), "\t", paste0(round(means,0), " ± ", round(stderrs,2), collapse="\t"))
}

getEffect(mcross, gtsp, 650, "all", "4_52450964")
getEffect(mcross, gtsp, 395, "all", "13_89689878")
getEffect(mcross, gtsp, 920, "all", "15_79013396")
getEffect(mcross, gtsp, 695, "all", "17_30971909")
getEffect(mcross, gtsp, 365, "all", "18_60822951")

getEffect(mcross, gtsp, 365, 1, "4_52450964", "longevity ~ site + cohort + treatment")
getEffect(mcross, gtsp, 755, 1, "6_140010438", "longevity ~ site + cohort + treatment")
getEffect(mcross, gtsp, 605, 1, "17_18001459", "longevity ~ site + cohort + treatment")






