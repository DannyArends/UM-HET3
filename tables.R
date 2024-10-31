setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
gtsp <- pull.geno(mcross)

pORm <- names(which(lapply(lapply(apply(gtsp,2,unique,na.rm = TRUE), na.omit), length) == 2))
gtsp <- gtsp[, pORm]

matM <- names(which(apply(gtsp,2, function(x){any(7 %in% x)})))
gtsp <- gtsp[, matM]
gtsp[gtsp == 7] <- "?C"
gtsp[gtsp == 8] <- "?D"


#### Paternal

lods.mM <- read.table("progressiveMapping_pat_males.txt", sep = "\t")
lods.fM <- read.table("progressiveMapping_pat_females.txt", sep = "\t")
lods.cM <- read.table("progressiveMapping_pat_all.txt", sep = "\t")


getEffect <- function(mcross, gtsprob, timepoint = 365, sex = "all", marker = "1_3010274", model = "longevity ~ sex + site + cohort + treatment"){
  mp <- gtsprob[, grep(marker, colnames(gtsprob))]

  cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]),
                      geno = mp)
  idx <- which(cdata[, "longevity"] >= timepoint)
  cat("N=", length(idx), "\n")
  cdata <- cdata[idx,]
  if(sex == 0 || sex == 1){ # males = 1, fem = 0
    cdata <- cdata[which(cdata[, "sex"] == sex),]
  }
  cat("N=", nrow(cdata), "\n")

  gts <- cdata[,"geno"]
  mlm <- lm(as.formula(model), data = cdata)

  pheAdj <- rep(NA, nrow(cdata))
  names(pheAdj) <-  rownames(cdata)

  adj <- residuals(mlm) + mean(cdata[, "longevity"])
  pheAdj[names(adj)] <- adj
  OAmean <- mean(pheAdj[which(!is.na(gts))])
  means <- c(mean(pheAdj[which(gts == "?C")]),mean(pheAdj[which(gts == "?D")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "?C")]),std(pheAdj[which(gts == "?D")]))
  paste0(length(which(!is.na(gts))), "/",nrow(cdata), "\t", round(OAmean,0), "\t", paste0(round(means,0), " ± ", round(stderrs,2), collapse="\t"))
}



getEffect(mcross, gtsp, marker = "1_132295971", timepoint = 365)

