#
# correlation.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Testing raw correlation between observed longevity and bodyweight measurements for combined, males, and females
#

library(svglite)

library(qtl)

source("ActuarialMapping/adjustXprobs.R")

mcross <- read.cross(format="csvr", file="DataSet/um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)
sex <- pull.pheno(mcross)[, "sex"]

# Include animals older than 365 days
above <- which(pull.pheno(mcross)[, "longevity"] >= 365)

observed <- pull.pheno(mcross)[above, "longevity"]
cor.test(observed, pull.pheno(mcross)[above, 5], use = "pair", method = "spearman")
cor.test(observed, pull.pheno(mcross)[above, 6], use = "pair", method = "spearman")
cor.test(observed, pull.pheno(mcross)[above, 7], use = "pair", method = "spearman")
cor.test(observed, pull.pheno(mcross)[above, 8], use = "pair", method = "spearman")

above <- which(pull.pheno(mcross)[, "longevity"] >= 365 & pull.pheno(mcross)[, "sex"] == 0)

observed <- pull.pheno(mcross)[above, "longevity"]
cor.test(observed, pull.pheno(mcross)[above, 5], use = "pair", method = "spearman")
cor.test(observed, pull.pheno(mcross)[above, 6], use = "pair", method = "spearman")
cor.test(observed, pull.pheno(mcross)[above, 7], use = "pair", method = "spearman")
cor.test(observed, pull.pheno(mcross)[above, 8], use = "pair", method = "spearman")

above <- which(pull.pheno(mcross)[, "longevity"] >= 365 & pull.pheno(mcross)[, "sex"] == 1)

observed <- pull.pheno(mcross)[above, "longevity"]
cor.test(observed, pull.pheno(mcross)[above, 5], use = "pair", method = "spearman")
cor.test(observed, pull.pheno(mcross)[above, 6], use = "pair", method = "spearman")
cor.test(observed, pull.pheno(mcross)[above, 7], use = "pair", method = "spearman")
cor.test(observed, pull.pheno(mcross)[above, 8], use = "pair", method = "spearman")


