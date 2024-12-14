setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

#
# Use a different clor scheme (blue - red)
#

### TODO: CTL for Males and Females separately (Done)

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

gtsp <- pull.genoprob(mcross)
cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    adjLongevity = NA, 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw24 = as.numeric(pull.pheno(mcross)[, "bw24"]),
                    adjBw24 = NA,
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
idx <- which(cdata[, "longevity"] >= 365 & !is.na(cdata[, "bw24"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw24 <- lm(bw24 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw24"] <- round(as.numeric(coef(lm.null.bw24)["(Intercept)"]) + residuals(lm.null.bw24), 2)

bCor <- cor(cdata[, "adjLongevity"], cdata[, "adjBw24"], use = "pair", method = "spearman")
bCor.f <- cor(cdata[which(cdata[, "sex"] == 0), "adjLongevity"], cdata[which(cdata[, "sex"] == 0), "adjBw24"], use = "pair", method = "spearman")
bCor.m <- cor(cdata[which(cdata[, "sex"] == 1), "adjLongevity"], cdata[which(cdata[, "sex"] == 1), "adjBw24"], use = "pair", method = "spearman")

