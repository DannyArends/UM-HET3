#
# Heritability_Main.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Compute Heritability using a full model of main Vita effects
# Uses an mean square methods (Adapted from: Falconer 1989 & Lynch & Walsh 1998)
#

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

markers <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085", "X_156343080")

names(markers) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa", "VitaXb")

rmarkers <- gtsp[,unlist(lapply(markers, grep, colnames(gtsp)))]

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]),
                    GenoID = as.character(pull.pheno(mcross)[, "GenoID"]),
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

sequ <-seq(35, 1100, 15)
reps <- 1:50

mmF <- vector("list", length(sequ))
names(mmF) <- sequ

tp <- 720
idx <- which(cdata[, "longevity"] >= tp)
adata <- cdata[idx, ]

Vita1a <- rmarkers[idx, 1:4]
Vita1b <- rmarkers[idx, 5:8]
Vita1c <- rmarkers[idx, 9:12]
Vita1d <- rmarkers[idx, 13:16]
Vita2a <- rmarkers[idx, 17:20]
Vita2b <- rmarkers[idx, 21:24]
Vita2c <- rmarkers[idx, 25:28]
Vita3a <- rmarkers[idx, 29:32]
Vita4a <- rmarkers[idx, 33:36]
Vita4b <- rmarkers[idx, 37:40]
Vita5a <- rmarkers[idx, 41:44]
Vita6a <- rmarkers[idx, 45:48]
Vita6b <- rmarkers[idx, 49:52]
Vita9a <- rmarkers[idx, 53:56]
Vita9b <- rmarkers[idx, 57:60]
Vita9c <- rmarkers[idx, 61:64]
Vita10a <- rmarkers[idx, 65:68]
Vita11a <- rmarkers[idx, 69:72]
Vita11b <- rmarkers[idx, 73:76]
Vita11c <- rmarkers[idx, 77:80]
Vita12a <- rmarkers[idx, 81:84]
Vita13a <- rmarkers[idx, 85:88]
Vita14a <- rmarkers[idx, 89:92]
Vita14b <- rmarkers[idx, 93:96]
Vita15a <- rmarkers[idx, 97:100]
Vita15b <- rmarkers[idx, 101:104]
Vita17a <- rmarkers[idx, 105:108]
Vita18a <- rmarkers[idx, 109:112]
VitaXa <- rmarkers[idx, 113:116]
VitaXb <- rmarkers[idx, 117:120]


mEff <- c()

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita1a = round(coef(lm(residuals(ml) ~ Vita1a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita1b = round(coef(lm(residuals(ml) ~ Vita1b + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita1c = round(coef(lm(residuals(ml) ~ Vita1c + 0)),1))


ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita1d = round(coef(lm(residuals(ml) ~ Vita1d + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita2a = round(coef(lm(residuals(ml) ~ Vita2a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2c + Vita3a + Vita4a + Vita4b + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita2b = round(coef(lm(residuals(ml) ~ Vita2b + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita3a + Vita4a + Vita4b + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita2c = round(coef(lm(residuals(ml) ~ Vita2c + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita4a + Vita4b + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita3a = round(coef(lm(residuals(ml) ~ Vita3a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4b + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita4a = round(coef(lm(residuals(ml) ~ Vita4a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita5a + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita4b = round(coef(lm(residuals(ml) ~ Vita4b + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita5a = round(coef(lm(residuals(ml) ~ Vita5a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita6a = round(coef(lm(residuals(ml) ~ Vita6a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita6b = round(coef(lm(residuals(ml) ~ Vita6b + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita9a = round(coef(lm(residuals(ml) ~ Vita9a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita9b = round(coef(lm(residuals(ml) ~ Vita9b + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita9c = round(coef(lm(residuals(ml) ~ Vita9c + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita10a = round(coef(lm(residuals(ml) ~ Vita10a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11b + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita11a = round(coef(lm(residuals(ml) ~ Vita11a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11c + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita11b = round(coef(lm(residuals(ml) ~ Vita11b + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita12a + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita11c = round(coef(lm(residuals(ml) ~ Vita11c + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita13a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita12a = round(coef(lm(residuals(ml) ~ Vita12a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + 
         Vita14a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita13a = round(coef(lm(residuals(ml) ~ Vita13a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + 
         Vita13a + Vita14b + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita14a = round(coef(lm(residuals(ml) ~ Vita14a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + 
         Vita13a + Vita14a + Vita15a  + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita14b = round(coef(lm(residuals(ml) ~ Vita14b + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + 
         Vita13a + Vita14a + Vita14b + Vita15b + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita15a = round(coef(lm(residuals(ml) ~ Vita15a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + 
         Vita13a + Vita14a + Vita14b + Vita15a + Vita17a + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita15b = round(coef(lm(residuals(ml) ~ Vita15b + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + 
         Vita13a + Vita14a + Vita14b + Vita15a + Vita15b + Vita18a + VitaXa, data = adata)

mEff <- rbind(mEff, vita17a = round(coef(lm(residuals(ml) ~ Vita17a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + 
         Vita13a + Vita14a + Vita14b + Vita15a + Vita15b + Vita17a + VitaXa, data = adata)

mEff <- rbind(mEff, vita18a = round(coef(lm(residuals(ml) ~ Vita18a + 0)),1))

ml <- lm(longevity ~ sex + site + cohort + treatment + 
         Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + 
         Vita5a + Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + 
         Vita13a + Vita14a + Vita14b + Vita15a + Vita15b + Vita17a + Vita18a, data = adata)

mEff <- rbind(mEff, vitaXa = round(coef(lm(residuals(ml) ~ VitaXa + 0)),1))


library(RColorBrewer)
colz <- brewer.pal(9, "PiYG")
colz <- colorRampPalette(colz)(80)

layout(matrix(c(1,1,1,1,1,1,1,1,2), ncol=9))

image(1:nrow(mEff), 1:ncol(mEff), mEff, xaxt="n", yaxt = "n", xlab = "", ylab = "", breaks = seq(-20, 20, 0.5), col = colz, 
  main = paste0("Individuals Vita effects at T", tp,", conditional on all other Vita loci"))
axis(1, at = 1:nrow(mEff), rownames(mEff), las = 2)
axis(2, at = 1:ncol(mEff), c("AC", "BC", "AD", "BD"), las=2)

image(1, seq(-20,20,0.5), t(seq(-20,20,0.5)), breaks = seq(-20, 20, 0.5), col = colz, xlab = "Effect (Days)", ylab = "",xaxt="n", las=2)
box()

