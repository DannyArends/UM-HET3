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

trait <- read.csv("ITP_10003.csv", header = TRUE, comment.char = "#", skip=10, row.names=2,na.strings = c("NA", "", "x"))
trait <- trait[which(trait[, "DA2024"] == 1),]
snames <- as.character(pull.pheno(mcross)[, "GenoID"])


gtsp <- pull.genoprob(mcross)
cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    adjLongevity = NA, 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw12 = as.numeric(trait[snames, "Value"]),
                    adjBw12 = NA,
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
rownames(cdata) <- pull.pheno(mcross)[, "GenoID"]

idx <- which(cdata[, "longevity"] >= 365 & !is.na(cdata[, "bw12"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw12 <- lm(bw12 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw12"] <- round(as.numeric(coef(lm.null.bw12)["(Intercept)"]) + residuals(lm.null.bw12), 2)

bCor <- cor(cdata[, "adjLongevity"], cdata[, "adjBw12"], use = "pair", method = "spearman")
bCor.f <- cor(cdata[which(cdata[, "sex"] == 0), "adjLongevity"], cdata[which(cdata[, "sex"] == 0), "adjBw12"], use = "pair", method = "spearman")
bCor.m <- cor(cdata[which(cdata[, "sex"] == 1), "adjLongevity"], cdata[which(cdata[, "sex"] == 1), "adjBw12"], use = "pair", method = "spearman")

