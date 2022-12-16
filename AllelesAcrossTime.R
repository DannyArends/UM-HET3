setwd("C:/Users/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Users/rqdt9/OneDrive - Northumbria University - Production Azure AD/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 1)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

# Our Progressive Mapping Sequence
msequence <- seq(365, 1100, 15)
markers <- unique(unlist(lapply(strsplit(colnames(gtsp), ":"), "[",1)))

lods.cM <- c()
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


#Combined
"1_3010272"                         , timepoint = 860)
"1_24042124"                         , timepoint = 695)
"2_112712327"                         , timepoint = 800)
"4_55012301"                         , timepoint = 650)
"6_107382038"                         , timepoint = 500)
"9_104091597"                         , timepoint = 1025)
"10_72780332"                         , timepoint = 980)
"12_112855820"                         , timepoint = 635)
"13_89689878"                         , timepoint = 395)
"14_101437457"                         , timepoint = 860)
"15_74248242"                         , timepoint = 905)
"17_32883804"                         , timepoint = 665)
"18_60822951"                         , timepoint = 365)
"X_36008085"                         , timepoint = 365)
"X_156343080"                         , timepoint = 740)

"1_24042124"                         , timepoint = 770, sex = 0, model = "longevity ~ sex + site + cohort")
"2_139956785"                         , timepoint = 545, sex = 0, model = "longevity ~ sex + site + cohort")
"3_92135706"                         , timepoint = 560, sex = 0, model = "longevity ~ sex + site + cohort")
"9_29939029"                         , timepoint = 665, sex = 0, model = "longevity ~ sex + site + cohort")
"9_104091597"                         , timepoint = 1040, sex = 0, model = "longevity ~ sex + site + cohort")
"11_82178599"                         , timepoint = 1040, sex = 0, model = "longevity ~ sex + site + cohort")
"X_156343080"                         , timepoint = 740, sex = 0, model = "longevity ~ sex + site + cohort")

"1_3010272"                         , timepoint = 890, sex = 1, model = "longevity ~ sex + site + cohort")
"1_120474787"                         , timepoint = 365, sex = 1, model = "longevity ~ sex + site + cohort")
"2_112712327"                         , timepoint = 935, sex = 1, model = "longevity ~ sex + site + cohort")
"3_83838529"                         , timepoint = 1070, sex = 1, model = "longevity ~ sex + site + cohort")
"4_52524395"                         , timepoint = 650, sex = 1, model = "longevity ~ sex + site + cohort")
"5_67573068"                         , timepoint = 1085, sex = 1, model = "longevity ~ sex + site + cohort")
"6_134870385"                         , timepoint = 365, sex = 1, model = "longevity ~ sex + site + cohort")
"9_124056586"                         , timepoint = 785, sex = 1, model = "longevity ~ sex + site + cohort")
"10_72780332"                         , timepoint = 980, sex = 1, model = "longevity ~ sex + site + cohort")
"11_5628810"                         , timepoint = 635, sex = 1, model = "longevity ~ sex + site + cohort")
"15_62405371"                         , timepoint = 845, sex = 1, model = "longevity ~ sex + site + cohort")
"17_34460077"                         , timepoint = 695, sex = 1, model = "longevity ~ sex + site + cohort")

