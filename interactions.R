#
# UM-HET3 interactions between significant Vita loci
#

all <- c("1_3010272", "1_24042124", "1_120474787", "2_112712327", "2_139956785", "2_157112564","3_83838529", "3_92135706", "4_55012301",
         "4_74811205", "4_145301445", "5_67573068", "6_54992703", "6_107382038", "8_36994142", "9_29939029", "9_54904313", "9_104091597",
         "9_124056586", "10_72780332", "11_5628810", "11_82178599", "12_112855820", "13_20905668", "13_89689878", "14_101437457", 
         "15_74248242", "17_32883804", "17_34460077", "17_68770703", "18_60822951", "X_36008085", "X_156343080")

names(all) <- c("Vita1a", "Vita1b", "Vita1c", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita3b", "Vita4a", "Vita4b", "Vita4c", "Vita5a",
                "Vita6a", "Vita6b", "Vita8a", "Vita9a", "Vita9b", "Vita9c", "Vita9d", "Vita10a", "Vita11a", "Vita11b", "Vita12a", "Vita13a",
                "Vita13b", "Vita14a", "Vita15a", "Vita17a", "Vita17b", "Vita17c", "Vita18a", "VitaXa", "VitaXb")

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

lodM <- c()
for(m1 in 1:length(all)){
  lodR <- c()
  for(m2 in 1:length(all)){
    mp1 <- gtsp[, grep(all[m1], colnames(gtsp))]
    mp2 <- gtsp[, grep(all[m2], colnames(gtsp))]
    lm.null <- lm(longevity ~ sex + site + cohort + treatment + mp1 + mp2 + 0, data = cdata)
    lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp1 * mp2 + 0, data = cdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lodR <- c(lodR, lod)
  }
  lodM <- rbind(lodM, lodR)
}
rownames(lodM) <- names(all)
colnames(lodM) <- names(all)

lodM[lower.tri(lodM)] <- NA
diag(lodM) <- NA

library(RColorBrewer)
colz <- c("white", brewer.pal(4,"PuRd"))

image(1:length(all), 1:length(all), lodM, xaxt="n", yaxt="n", xlab="", ylab="", breaks = c(0, 3, 6, 9, 12, 100), col=colz)
axis(1, at = 1:length(all), names(all), las=2)
axis(2, at = 1:length(all), names(all), las=2)

write.table(lodM, "vita_interactions.txt", sep = "\t", quote=FALSE)


