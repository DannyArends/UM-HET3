setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

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
msequence <- seq(20, 1100, 60)
markers <- unique(unlist(lapply(strsplit(colnames(gtsp), ":"), "[",1)))

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny/GxD")

lods.cM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.factor(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] >= x)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.null <- lm(longevity ~ sex + site + cohort + treatment + mp + 0, data = cdata)
    lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp + mp:treatment + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.cM <- rbind(lods.cM, lods.c)
  cat("Done", x, "\n")
}
colnames(lods.cM) <- colnames(pull.geno(mcross))
rownames(lods.cM) <- paste0("> ", msequence)
write.table(lods.cM, "combined_GxD.txt", sep = "\t", quote = FALSE)

lods.fM <- c()
minAge <- c()
for(x in msequence){
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.factor(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] >= x & cdata[, "sex"] == 0)
  cdata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  minAge <- c(minAge, min(cdata[, "longevity"]))

  lods.c <- c()
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.null <- lm(longevity ~ site + cohort + treatment + mp + 0, data = cdata)
    lm.alt <- lm(longevity ~ site + cohort + treatment + mp + mp:treatment + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.fM <- rbind(lods.fM, lods.c)
  cat("Done", x, "\n")
}
colnames(lods.fM) <- colnames(pull.geno(mcross))
rownames(lods.fM) <- paste0("> ", msequence)
write.table(lods.fM, "female_GxD.txt", sep = "\t", quote = FALSE)


lods.mM <- c()
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
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.null <- lm(longevity ~ site + cohort + treatment + mp + 0, data = cdata)
    lm.alt <- lm(longevity ~ site + cohort + treatment + mp + mp:treatment + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.mM <- rbind(lods.mM, lods.c)
  cat("Done", x, "\n")
}
colnames(lods.mM) <- colnames(pull.geno(mcross))
rownames(lods.mM) <- paste0("> ", msequence)
write.table(lods.mM, "male_GxD.txt", sep = "\t", quote = FALSE)

