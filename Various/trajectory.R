#
# trajectory.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Bodyweight trajectory analysis (DeCabo)
#

# Known Loci
allV <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085")

names(allV) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa")

allS <- c("1_3010272","1_86216552","2_13600088", "2_60201233","2_161871392","3_87974845","3_159581164","4_30761996","4_107374161","6_8006720",
         "6_138658041","7_16072018","7_120086292","8_71684276","8_126505019","9_51116640","10_18144599","11_97448477","12_71677220",
         "12_118179607","13_19367506","13_98521647","14_30957748", "14_101437466","15_3288506","16_75758401","17_26542857",
         "18_52488251","19_3403302","19_53851357")
names(allS) <- c("Soma1a","Soma1b","Soma2a","Soma2b","Soma2c","Soma3a","Soma3b","Soma4a","Soma4b","Soma6a","Soma6b","Soma7a","Soma7b","Soma8a",
                "Soma8b", "Soma9a", "Soma10a","Soma11a","Soma12a","Soma12b","Soma13a","Soma13b","Soma14a","Soma14b","Soma15a","Soma16a",
                "Soma17a","Soma18a","Soma19a","Soma19b")

all <- c(allV, allS)

setwd("C:/Github/UM-HET3/ProgessiveMapping")
source("adjustXprobs.R")
setwd("C:/Users/Danny/OneDrive - Northumbria University - Production Azure AD/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)
sex <- pull.pheno(mcross)[, "sex"]
rownames(gtsp) <- pull.pheno(mcross)[,"GenoID"]

# Read Census & Subset
setwd("C:/Users/Danny/UTHSC GGI Dropbox/Danny Arends/MyFolder/UM-HET3/DeCabo")
library(nnet)
lodM <- c()

computeLog10 <- function(lm.null, lm.alt){
  # Extract Log-Likelihood values
  LL_null <- logLik(lm.null)[1]
  LL_alt <- logLik(lm.alt)[1]
  G2_manual <- 2 * (LL_alt - LL_null)
  df_null <- length(coef(lm.null))
  df_alt <- length(coef(lm.alt))
  df_manual <- df_alt - df_null

  # Compute the P-value using the Chi-Squared distribution
  return(-log10(pchisq(G2_manual, df = df_manual, lower.tail = FALSE)))
}


### File 1 - Census
census <- read.csv("itp_geno_census.csv", row.names = 15)
gtsM <- gtsp[rownames(census),]

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.factor(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
rownames(cdata) <- pull.pheno(mcross)[,"GenoID"]
cdata <- cdata[rownames(census),]
cdata <- cbind(new_class_bw = as.factor(census[, "new_class_bw"]), cdata)

lods.c <- c()
lods.mn <- c()
lm.null <- multinom(new_class_bw ~ sex + site + cohort + 0, data = cdata, trace = FALSE)
for(marker in all){
  mp <- gtsM[, grep(marker, colnames(gtsM))]
  lm.alt <- multinom(new_class_bw ~ sex + site + cohort + mp + 0, data = cdata, trace = FALSE)
  lods.c <- c(lods.c, computeLog10(lm.null, lm.alt))
}
names(lods.c) <- names(all)
lodM <- cbind(lodM, lods.c)

### File 2 - Census F
census <- read.csv("itp_geno_f_census.csv", row.names = 11)
gtsM <- gtsp[rownames(census),]

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.factor(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
rownames(cdata) <- pull.pheno(mcross)[,"GenoID"]
cdata <- cdata[rownames(census),]
cdata <- cbind(new_class_bw = census[, "new_class_bw"], cdata)

lods.c <- c()
lm.null <- multinom(new_class_bw ~ site + cohort + 0, data = cdata, trace = FALSE)
for(marker in all){
  mp <- gtsM[, grep(marker, colnames(gtsM))]
  lm.alt <- multinom(new_class_bw ~ site + cohort + mp + 0, data = cdata, trace = FALSE)
  lods.c <- c(lods.c, computeLog10(lm.null, lm.alt))
}
names(lods.c) <- names(all)
lodM <- cbind(lodM, lods.c)

### File 3 - Census M
census <- read.csv("itp_geno_m_census.csv", row.names = 12)
gtsM <- gtsp[rownames(census),]

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.factor(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
rownames(cdata) <- pull.pheno(mcross)[,"GenoID"]
cdata <- cdata[rownames(census),]
cdata <- cbind(new_class_bw = census[, "new_class_bw"], cdata)

lods.c <- c()
lm.null <- multinom(new_class_bw ~ site + cohort + 0, data = cdata, trace = FALSE)
for(marker in all){
  mp <- gtsM[, grep(marker, colnames(gtsM))]
  lm.alt <- multinom(new_class_bw ~ site + cohort + mp + 0, data = cdata, trace = FALSE)
  lods.c <- c(lods.c, computeLog10(lm.null, lm.alt))
}
names(lods.c) <- names(all)
lodM <- cbind(lodM, lods.c)


### File 4 - Census TX
census <- read.csv("itp_geno_tx_census.csv", row.names = 17)
gtsM <- gtsp[rownames(census),]

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.factor(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
rownames(cdata) <- pull.pheno(mcross)[,"GenoID"]
cdata <- cdata[rownames(census),]
cdata <- cbind(new_class_bw = census[, "new_class_bw"], cdata)

lods.c <- c()
lm.null <- multinom(new_class_bw ~ sex + site + cohort + treatment + 0, data = cdata, trace = FALSE)
for(marker in all){
  mp <- gtsM[, grep(marker, colnames(gtsM))]
  lm.alt <- multinom(new_class_bw ~ sex + site + cohort + treatment + mp + 0, data = cdata, trace = FALSE)
  lods.c <- c(lods.c, computeLog10(lm.null, lm.alt))
}
names(lods.c) <- names(all)
lodM <- cbind(lodM, lods.c)

### File 5 - Census TX F
census <- read.csv("itp_geno_tx_f_census.csv", row.names = 13)
gtsM <- gtsp[rownames(census),]

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.factor(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
rownames(cdata) <- pull.pheno(mcross)[,"GenoID"]
cdata <- cdata[rownames(census),]
cdata <- cbind(new_class_bw = census[, "new_class_bw"], cdata)

lods.c <- c()
lm.null <- multinom(new_class_bw ~ site + cohort + treatment + 0, data = cdata, trace = FALSE)
for(marker in all){
  mp <- gtsM[, grep(marker, colnames(gtsM))]
  lm.alt <- multinom(new_class_bw ~ site + cohort + treatment + mp + 0, data = cdata, trace = FALSE)
  lods.c <- c(lods.c, computeLog10(lm.null, lm.alt))
}
names(lods.c) <- names(all)
lodM <- cbind(lodM, lods.c)


### File 6 - Census TX M
census <- read.csv("itp_geno_tx_m_census.csv", row.names = 13)
gtsM <- gtsp[rownames(census),]

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.factor(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
rownames(cdata) <- pull.pheno(mcross)[,"GenoID"]
cdata <- cdata[rownames(census),]
cdata <- cbind(new_class_bw = census[, "new_class_bw"], cdata)

lods.c <- c()
lm.null <- multinom(new_class_bw ~ site + cohort + treatment + 0, data = cdata, trace = FALSE)
for(marker in all){
  mp <- gtsM[, grep(marker, colnames(gtsM))]
  lm.alt <- multinom(new_class_bw ~ site + cohort + treatment + mp + 0, data = cdata, trace = FALSE)
  lods.c <- c(lods.c, computeLog10(lm.null, lm.alt))
}
names(lods.c) <- names(all)
lodM <- cbind(lodM, lods.c)

lodM <- round(lodM, 2)
colnames(lodM) <- c("census", "census_f", "census_m", "census_tx", "census_tx_f", "census_tx_m")

chr <- as.character(lapply(strsplit(all,"_"),"[",1))
pos <- as.character(lapply(strsplit(all,"_"),"[",2))
lodA <- cbind(chr, pos, lodM)
write.table(lodA, "census_mapping.txt", sep = "\t", quote = FALSE)