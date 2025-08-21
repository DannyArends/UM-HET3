setwd("D:/Ddrive/Github/UM-HET3")
source("adjustXprobs.R")

# Change setwd() to where you stored your data
setwd("C:/Github/UM-HET3/files")

library(qtl)

# Read cross object & compute genoprob
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

### Quick n dirty R/QTL scan, using sex, site, cohort and treatment as covariates
covariates <- data.frame(sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                         site = as.factor(pull.pheno(mcross)[, "site"]), 
                         treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

addcovar <- model.matrix(~ sex + site + treatment + 0, covariates)

allI <- scanone(mcross, pheno.col = "longevity", addcovar=addcovar[,-2], method="hk")

setwd("D:/Ddrive/GN2fixCovar")
zcross <- read.cross(format="csvr", file="4B8Oj0tF.cross", genotypes=NULL, na.strings=c("-", "NA"))
zcross <- calc.genoprob(zcross)
zcross <- adjustXprobs(zcross)

### Quick n dirty R/QTL scan, using sex, site, cohort and treatment as covariates
covariates <- data.frame(sex = as.numeric(pull.pheno(zcross)[, "T_10009"]), 
                         site = as.factor(pull.pheno(zcross)[, "T_10011"]), 
                         treatment = as.factor(pull.pheno(zcross)[, "T_10012"]))

addcovar <- model.matrix(~ sex + site + treatment + 0, covariates)

allZ <- scanone(zcross, pheno.col = "T_10001", addcovar=addcovar[,-2], method="hk")



### Older TESTS

### R/QTL scan using 1 covariate (sex)
covariates <- data.frame(sex = as.numeric(pull.pheno(mcross)[, "sex"]))
res.rqtl <- scanone(mcross, addcovar = covariates, pheno.col = "longevity", method="hk")

### Use the R/qtl computed genotype probabilities and do linear regression using these probabilities
gtsp <- pull.genoprob(mcross)

longevity <- pull.pheno(mcross)[, "longevity"]
sex <- as.numeric(pull.pheno(mcross)[, "sex"])

lods <- c()
for(marker in colnames(pull.geno(mcross))){
  mp <- gtsp[,grep(marker, colnames(gtsp))]
  lod <- -log10(anova(lm(longevity ~ sex + mp + 1), lm(longevity ~ sex + 1))[, "Pr(>F)"][2])
  lods <- c(lods, lod)
}
res.gplm <- res.rqtl
res.gplm[, 3] <- NA
res.gplm[colnames(pull.geno(mcross)), 3] <- lods

# Compare, everything look similar with only minor differences in LOD
plot(res.rqtl, res.gplm, main="Comparison", col = c("black", "green"), lty=c(1, 2))
legend("topright", c("R/qtl", "LM using GT prob"), col = c("black", "green"), lty=c(1, 2))

######## Now using 2 covariates, 1 numeric & 1 factor

### R/QTL scan
covariates <- data.frame(sex = as.numeric(pull.pheno(mcross)[, "sex"]),
                         site = as.factor(pull.pheno(mcross)[, "site"]))

#Use the model.matrix function to 'explode' site into 3 columns
addcovar <- model.matrix(~ sex + site + 0, covariates)

res.rqtl <- scanone(mcross, addcovar = addcovar[,-2], pheno.col = "longevity", method="hk")

### Use the R/qtl computed genotype probabilities and do linear regression using these probabilities
gtsp <- pull.genoprob(mcross)

longevity <- pull.pheno(mcross)[, "longevity"]
sex <- as.numeric(pull.pheno(mcross)[, "sex"])
site <- as.factor(pull.pheno(mcross)[, "site"])
lods <- c()

for(marker in colnames(pull.geno(mcross))){
  mp <- gtsp[,grep(marker, colnames(gtsp))]
  lod <- -log10(anova(lm(longevity ~ sex + site + mp + 1), lm(longevity ~ sex + site + 1))[, "Pr(>F)"][2])
  lods <- c(lods, lod)
}
res.gplm <- res.rqtl
res.gplm[, 3] <- NA
res.gplm[colnames(pull.geno(mcross)), 3] <- lods

plot(res.rqtl, res.gplm, main="Comparison", col = c("black", "blue"), lwd=c(1,2), lty=c(1, 2))
legend("topright", c("R/qtl", "LM using GT prob"), col = c("black", "blue"), lwd=c(1,2), lty=c(1, 2))




### Quick n dirty R/QTL scan, using sex, site, cohort and treatment as covariates
covariates <- data.frame(sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                         site = as.factor(pull.pheno(mcross)[, "site"]), 
                         cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                         treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

addcovar <- model.matrix(~ sex + site + treatment + 0, covariates)

allI <- scanone(mcross, pheno.col = "longevity", addcovar=addcovar[,-2], method="hk")

plot(allI, main = "Longevity")
abline(h = -log10(0.01 / 400), col="Green")
abline(h = -log10(0.05 / 400), col="Orange")
abline(h = -log10(0.1 / 400), col="Red")
legend("topright",legend = c(0.1, 0.05, 0.01), col = c("red", "orange", "green"), lty=1)

# Subset males and females
cross.f <- subset(mcross, ind = pull.pheno(mcross)[, "sex"] == 0)
cross.m <- subset(mcross, ind = pull.pheno(mcross)[, "sex"] == 1)

covar.f <- cbind(as.numeric(pull.pheno(cross.f)[, "site"]), as.numeric(pull.pheno(cross.f)[, "cohort"], pull.pheno(cross.f)[, "treatment"]))
scan.f <- scanone(cross.f, pheno.col = "longevity", addcovar=covar.f, method="hk")

covar.m <- cbind(as.numeric(pull.pheno(cross.m)[, "site"]), as.numeric(pull.pheno(cross.m)[, "cohort"], pull.pheno(cross.m)[, "treatment"]))
scan.m <- scanone(cross.m, pheno.col = "longevity", addcovar=covar.m, method="hk")

plot(allI, scan.f, scan.m, main = "Longevity", col=c(1,"pink","blue"))
abline(h = -log10(0.01 / 400), col="Green")
abline(h = -log10(0.05 / 400), col="Orange")
abline(h = -log10(0.1 / 400), col="Red")
legend("topleft", legend=c("All", "Female", "Male"), lwd=2, col=c(1,"pink","blue"))

### Mapping using LMs

bw6m <- indM[, "bw6"]
bw12m <- indM[, "bw12"]
bw18m <- indM[, "bw18"]
bw24m <- indM[, "bw24"]
longevity <- indM[, "longevity"]

sex <- as.factor(indM[, "sex"])
site <- as.factor(indM[, "site"])
year <- as.factor(indM[, "cohort"])
treatment <- as.factor(indM[, "treatment"])

anova(lm(bw6m ~ sex + site + year + treatment)) # Sex, Site, Year
anova(lm(bw12m ~ sex + site + year + treatment)) # Sex, Site, Year, Treatment
anova(lm(bw18m ~ sex + site + year + treatment)) # Sex, Site, Year, Treatment
anova(lm(bw24m ~ sex + site + year + treatment)) # Sex, Site, Year
m1 <- lm(longevity ~ sex + site + year + treatment)
m2 <- lm(longevity ~ sex + site + year + treatment + sex:site + site:treatment)
AIC(m1,m2)
mAnova <- anova(lm(longevity ~ sex * site * year * treatment))
