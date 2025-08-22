#
# deathsPer15d.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Compute deaths in each 15 day window per genotype (>85% certain) per sex
# Additionally, at the end we create a visualization of the data
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

all <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085")

names(all) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa")

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__bioRxiv_All_Key_Files/11_FiguresDanny/")
for(pname in names(all)){
  marker <-  all[pname]

  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  longevity = as.numeric(pull.pheno(mcross)[, "longevity"])

  ii1 <- which(gts == "AC") # CH
  ii2 <- which(gts == "AD") # CD
  ii3 <- which(gts == "BC") # BH
  ii4 <- which(gts == "BD") # BD

  nnCH <- length(ii1)
  nnCD <- length(ii2)
  nnBH <- length(ii3)
  nnBD <- length(ii4)

  mm <- c()
  for(x in seq(20, 1500, 15)){
    iiD <- which(longevity >= x-15 & longevity < x)
    nCH <- length(which(iiD %in% ii1))
    nCD <- length(which(iiD %in% ii2))
    nBH <- length(which(iiD %in% ii3))
    nBD <- length(which(iiD %in% ii4))
    mm <- rbind(mm, c(nCH, nCD, nBH, nBD, nCH + nCD, nBH + nBD, nCD + nBD, nCH + nBH))
  }
  colnames(mm) <- c("CH", "CD", "BH", "BD", "C", "B", "D", "H")
  rownames(mm) <- paste0(seq(5, 1500, 15), "-", seq(20, 1515, 15))[-100]
  write.table(mm, paste0(pname, "_deathsIn15Dwindows_combined.txt"), sep = "\t", quote = FALSE)
}

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__bioRxiv_All_Key_Files/11_FiguresDanny/")
for(pname in names(all)){
  marker <-  all[pname]

  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  longevity = as.numeric(pull.pheno(mcross)[, "longevity"])
  sex = as.numeric(pull.pheno(mcross)[, "sex"])

  ii1 <- which(gts == "AC" & sex == 0) # CH
  ii2 <- which(gts == "AD" & sex == 0) # CD
  ii3 <- which(gts == "BC" & sex == 0) # BH
  ii4 <- which(gts == "BD" & sex == 0) # BD

  nnCH <- length(ii1)
  nnCD <- length(ii2)
  nnBH <- length(ii3)
  nnBD <- length(ii4)

  mm <- c()
  for(x in seq(20, 1500, 15)){
    iiD <- which(longevity >= x-15 & longevity < x & sex == 0)
    nCH <- length(which(iiD %in% ii1))
    nCD <- length(which(iiD %in% ii2))
    nBH <- length(which(iiD %in% ii3))
    nBD <- length(which(iiD %in% ii4))
    mm <- rbind(mm, c(nCH, nCD, nBH, nBD, nCH + nCD, nBH + nBD, nCD + nBD, nCH + nBH))
  }
  colnames(mm) <- c("CH", "CD", "BH", "BD", "C", "B", "D", "H")
  rownames(mm) <- paste0(seq(5, 1500, 15), "-", seq(20, 1515, 15))[-100]
  write.table(mm, paste0(pname, "_deathsIn15Dwindows_female.txt"), sep = "\t", quote = FALSE)
}

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__bioRxiv_All_Key_Files/11_FiguresDanny/")
for(pname in names(all)){
  marker <-  all[pname]

  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  longevity = as.numeric(pull.pheno(mcross)[, "longevity"])
  sex = as.numeric(pull.pheno(mcross)[, "sex"])

  ii1 <- which(gts == "AC" & sex == 1) # CH
  ii2 <- which(gts == "AD" & sex == 1) # CD
  ii3 <- which(gts == "BC" & sex == 1) # BH
  ii4 <- which(gts == "BD" & sex == 1) # BD

  nnCH <- length(ii1)
  nnCD <- length(ii2)
  nnBH <- length(ii3)
  nnBD <- length(ii4)

  mm <- c()
  for(x in seq(20, 1500, 15)){
    iiD <- which(longevity >= x-15 & longevity < x & sex == 1)
    nCH <- length(which(iiD %in% ii1))
    nCD <- length(which(iiD %in% ii2))
    nBH <- length(which(iiD %in% ii3))
    nBD <- length(which(iiD %in% ii4))
    mm <- rbind(mm, c(nCH, nCD, nBH, nBD, nCH + nCD, nBH + nBD, nCD + nBD, nCH + nBH))
  }
  colnames(mm) <- c("CH", "CD", "BH", "BD", "C", "B", "D", "H")
  rownames(mm) <- paste0(seq(5, 1500, 15), "-", seq(20, 1515, 15))[-100]
  write.table(mm, paste0(pname, "_deathsIn15Dwindows_male.txt"), sep = "\t", quote = FALSE)
}








#
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny")
#pdf("KM-Vita1A-DeathP15DWindow.pdf")
#par(cex=1.5)
#par(cex.axis=1.2)
plot(c(0, 1200), c(0.1, 50), t = "n", main = "Vita1a -Combined - % Deaths per Window", 
  ylab = "% of animals dead", xlab = "Age", yaxt = "n", log = "y")
mm <- c()

X <- seq(15, 1200, 15)

for(x in X){
  iiD <- which(longevity >= x-15 & longevity < x)
  nCH <- length(which(iiD %in% ii1))
  nBH <- length(which(iiD %in% ii2))
  nCD <- length(which(iiD %in% ii3))
  nBD <- length(which(iiD %in% ii4))
  cat(x, "", nCH, nBH, nCD, nBD, "\n")
  mm <- rbind(mm, c(nCH, nBH, nCD, nBD))
  #points(x, nAC, col = 1)
  #points(x, nAD, col = 2)
  #points(x, nBC, col = 3)
  #points(x, nBD, col = 4)
}
colnames(mm) <- c("CH", "BH", "CD", "BD")
rownames(mm) <- X
write.table(mm, "deathsIn15Dwindows.txt", sep = "\t", quote = FALSE)
Y1 <- 100 *mm[,1] / (nnCH - cumsum(mm[,1]))
#points(X, predict(loess(Y1 ~ X + 0)), t = "l", col = col.main[1])
points(X, Y1, t = "p", col = col.main[1])
abline(lm(log(Y1) ~ X), col = col.main[1])
Y2 <- 100 *mm[,2]/ (nnBH - cumsum(mm[,2]))
points(X, Y2, t = "p", col = col.main[2])
#points(X, predict(loess(Y2 ~ X + 0)), t = "l", col = col.main[2])
Y3 <- 100 *mm[,3]/ (nnCD- cumsum(mm[,3]))
points(X, Y3, t = "p", col = col.main[3])
#points(X, predict(loess(Y3 ~ X + 0)), t = "l", col = col.main[3])
Y4 <- 100 *mm[,4]/ (nnBD- cumsum(mm[,4]))
points(X, Y4, t = "p", col = col.main[4])
#points(X, predict(loess(Y4 ~ X + 0)), t = "l", col = col.main[4])
#dev.off()
axis(2, at = seq(0, 50, 2.5), seq(0, 50, 2.5),las=2)



nnCH <- length(ii1)
nnBH <- length(ii2)
nnCD <- length(ii3)
nnBD <- length(ii4)


pCH <- nnCH / 20 
pBH <- nnBH / 20 
pCD <- nnCD / 20 
pBD <- nnBD / 20 

round(seq(0, nnCH, pCH))
#setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny")
#pdf("TimeFor5percToDie.pdf")
#par(cex=1.5)
#par(cex.axis=1.2)
plot(c(1, 20), c(10, 500), t = "n", log = "y", main = "Time for 5% of the population to Die", 
     xaxt = "n", yaxt = "n", xlab = "5% mortality window", ylab = "Mortality window duration (d)")
points(diff(c(0, sort(longevity[ii1])[round(seq(0, nnCH, pCH))])), t = "l", col = col.main[1])
points(diff(c(0, sort(longevity[ii2])[round(seq(0, nnBH, pBH))])), t = "l", col = col.main[2])
points(diff(c(0, sort(longevity[ii3])[round(seq(0, nnCD, pCD))])), t = "l", col = col.main[3])
points(diff(c(0, sort(longevity[ii4])[round(seq(0, nnBD, pBD))])), t = "l", col = col.main[4])
axis(1, at = seq(1,20,1), seq(5, 100, 5))
axis(2, at = seq(0,500,10), seq(0, 500, 10), las=2)

#dev.off()
