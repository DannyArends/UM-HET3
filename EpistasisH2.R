setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]),
                    GenoID = as.character(pull.pheno(mcross)[, "GenoID"]),
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

sequ <-seq(35, 1100, 15)
reps <- 1:50

all <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085")

names(all) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", 
                "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", "Vita9a", "Vita9b", "Vita9c", 
                "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", 
                "Vita14b", "Vita15a", "Vita15b", "Vita17a", "Vita18a", "VitaXa")

rmarkers <- gtsp[,unlist(lapply(all, grep, colnames(gtsp)))]

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__bioRxiv_All_Key_Files/11_FiguresDanny/GxG_April25")
lodM.f <- read.table(paste0("vita_interactions_2way_females_tp42.txt"), sep = "\t")
lodM.m <- read.table(paste0("vita_interactions_2way_males_tp42.txt"), sep = "\t")

ints <- c()
for(timepoint in c(42, 365, 740, 905)){
  for(mX1 in 1:(length(all)-1)){
    for(mX2 in mX1:length(all)){
      m1 <- names(all)[mX1]
      m2 <- names(all)[mX2]
      if(m1 == m2) next;

      if(lodM.m[m2,m1] >= 3.8 || lodM.f[m2,m1] >= 3.8){
       ints <- rbind(ints, c(m1,m2, lodM.m[m2,m1] >= 3.8, lodM.f[m2,m1] >= 3.8))
       
      }
}}}
ints <- unique(ints)
colnames(ints) <- c("M1", "M2", "Males", "Females")
isMale <- ints[which(ints[,3] == "TRUE"),]
isFem <- ints[which(ints[,4] == "TRUE"),]
for(x in 1:nrow(isFem)){
  cat(paste0(isFem[x,1], ":", isFem[x,2]), " + ")
}

for(x in 1:nrow(isMale)){
  cat(paste0(isMale[x,1], ":", isMale[x,2]), " + ")
}
ints


mmF <- vector("list", length(sequ))
names(mmF) <- sequ
for(x in sequ){
  Hs <- c()
  set.seed(42)
  for(n in reps) {
    idx <- which(cdata[, "longevity"] >= x & cdata[, "sex"] == 0)
    idx <- sample(idx, length(idx) * .9)
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
    ml <- lm(longevity ~ site + cohort + treatment + 
             Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + Vita5a + 
             Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
             Vita14a + Vita14b + Vita15a + Vita15b  + Vita17a + Vita18a + VitaXa +
             Vita1a:Vita6b + Vita1a:Vita9a + Vita1a:Vita11a + Vita1b:Vita6b + Vita1b:Vita9a + 
             Vita1b:Vita9b + Vita1c:Vita2b + Vita1d:Vita2b + Vita1d:Vita15b + Vita1d:Vita18a + 
             Vita2b:Vita6a + Vita3a:Vita5a + Vita4b:Vita12a + Vita5a:Vita17a + Vita6a:Vita11b + 
             Vita6b:Vita14a + Vita6b:Vita14b + Vita11a:Vita14a + Vita11b:Vita14a
    , data = adata)
    ma <- anova(ml)
    # Mean square methods (most common one used)
    N <- (length(idx)/4)
    sQEE <- (ma[, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N            # Partial Total variance
    sQTT <- (ma["Residuals", "Mean Sq"] + sum(ma[-nrow(ma), "Mean Sq"])) / N    # Fixed effects & residual variance
    vE <- sQEE / sQTT # Contributions
    Hs <- cbind(Hs, vE)
  }
  cat("Done ", x,"\n")
  colnames(Hs) <- reps
  rownames(Hs) <- c("site", "cohort", "treatment", c(names(all), paste0("I", 1:19)), "unexplained")

  mmF[[as.character(x)]] <- Hs
}


mG <- c()
mE <- c()
for(x in sequ[-length(sequ)]){
  aa <- round(apply(mmF[[as.character(x)]][4:51,],1, mean),3)
  mG <- rbind(mG, aa)
  aa <- round(apply(mmF[[as.character(x)]][1:3,], 1, mean),3)
  mE <- rbind(mE, aa)
}
mG[mG < 0] <- 0
mE[mE < 0] <- 0
rownames(mG) <- sequ[-length(sequ)]
colnames(mG) <- c(names(all), paste0("I", 1:19))

rownames(mE) <- sequ[-length(sequ)]
colnames(mE) <- c("Site", "Cohort", "Treatment")


add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__bioRxiv_All_Key_Files/11_FiguresDanny")
pdf("H2_G_GxG_Female.pdf", width = 26, height = 26/3)
op <- par(mar = c(4,10,2,2))
op <- par(cex = 1.5)

plot(c(10, 1110), y = c(0, 0.50), xlab = "Tage (Days)", ylab = "Heritability", 
     main = "Female Heritability - 29 Vita regions + 19 GxG interactions", t = "n", yaxt="n", xaxt="n", yaxs = "i", xaxs = "i")
unexp <- c()
for(x in 1:(length(sequ)-1)){
  unexp <- c(unexp, 1 - (sum(mG[x,]) + sum(mE[x,])))
}

for(x in sequ){
  aa <- mmF[[as.character(x)]][4:32,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#FF3333", 0.4), medcol="#FF3333", border="#FF3333", outpch = 19, outcex = .5)

  aa <- mmF[[as.character(x)]][33:51,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#3333FF", 0.4), medcol="#3333FF", border="#3333FF", outpch = 19, outcex = .5)


  aa <- mmF[[as.character(x)]][1:3,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#466D1D", 0.4), medcol="#607D3B", border="#607D3B", outpch = 19, outcex = .5)
}

lo <- loess(unexp ~ sequ[-length(sequ)])
points(sequ[-length(sequ)], predict(lo), t = "l", col = "orange", lwd=3)
axis(1, at = seq(20, 1100, 4*15), seq(20, 1100, 4*15))
axis(2, at = seq(0, 1,0.1), paste0(seq(0,100,10), "%"), las=2)
axis(2, at = seq(0, 1,0.05), rep("", length(seq(0, 1,0.05))), las=2)
#legend("topleft", c("Female (G)", "Environment", "Female (GxG)"), fill = add.alpha(c("#FF3333", "#AFE1AF", "#3333FF"), 0.4), border = c("#FF3333", "#AFE1AF", "#3333FF"))
dev.off()


mmM <- vector("list", length(sequ))
names(mmM) <- sequ
for(x in sequ){
  Hs <- c()
  set.seed(42)
  for(n in reps){
    idx <- which(cdata[, "longevity"] >= x & cdata[, "sex"] == 1)
    idx <- sample(idx, length(idx) * .9)
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
    ml <- lm(longevity ~ site + cohort + treatment + 
             Vita1a + Vita1b + Vita1c + Vita1d + Vita2a + Vita2b + Vita2c + Vita3a + Vita4a + Vita4b + Vita5a + 
             Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita11c + Vita12a + Vita13a + 
             Vita14a + Vita14b + Vita15a + Vita15b  + Vita17a + Vita18a + VitaXa +
             Vita1a:Vita11c + Vita1c:Vita3a + Vita1c:Vita11c + Vita1d:Vita2c + Vita1d:Vita9b + 
             Vita1d:Vita9c + Vita2a:Vita13a + Vita2a:Vita14a + Vita2b:Vita3a + Vita2b:Vita14a + 
             Vita2b:Vita14b + Vita2c:Vita9b + Vita4b:Vita14b + Vita6a:Vita14a + Vita6b:Vita13a + 
             Vita6b:Vita15b + Vita9b:Vita17a + Vita9c:Vita17a + Vita10a:Vita14a + Vita11b:Vita17a + 
             Vita11c:Vita13a + Vita14a:Vita17a, data = adata);
    ma <- anova(ml)
    # Mean square methods (most common one used)
    N <- (length(idx)/4)
    sQEE <- (ma[, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N # Partial Total variance
    sQTT <- (ma["Residuals", "Mean Sq"] + sum(ma[-nrow(ma), "Mean Sq"])) / N # Fixed effects & residual variance
    vE <- sQEE / sQTT # Contributions
    Hs <- cbind(Hs, vE)
  }
  colnames(Hs) <- reps
  rownames(Hs) <- c("site", "cohort", "treatment", c(names(all), paste0("I", 1:22)), "unexplained")
  cat("Done ", x,"\n")
  mmM[[as.character(x)]] <- Hs
}



mG <- c()
mE <- c()
for(x in sequ){
  aa <- round(apply(mmM[[as.character(x)]][4:54,],1, mean),3)
  mG <- rbind(mG, aa)
  aa <- round(apply(mmM[[as.character(x)]][1:3,],1, mean),3)
  mE <- rbind(mE, aa)
}
mG[mG < 0] <- 0
mE[mE < 0] <- 0
rownames(mG) <- sequ[-length(sequ)]
colnames(mG) <- c(names(all), paste0("I", 1:22))

rownames(mE) <- sequ[-length(sequ)]
colnames(mE) <- c("Site", "Cohort", "Treatment")



add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__bioRxiv_All_Key_Files/11_FiguresDanny")
pdf("H2_G_GxG_Male.pdf", width = 26, height = 26/3)
op <- par(mar = c(4,10,2,2))
op <- par(cex = 1.5)

plot(c(10, 1110), y = c(0, 0.50), xlab = "Tage (Days)", ylab = "Heritability", 
     main = "Male Heritability - 29 Vita regions + 22 GxG interactions", t = "n", yaxt="n", xaxt="n", yaxs = "i", xaxs = "i")
unexp <- c()
for(x in 1:(length(sequ)-1)){
  unexp <- c(unexp, 1 - (sum(mG[x,]) + sum(mE[x,])))
}

for(x in sequ){
  aa <- mmM[[as.character(x)]][4:32,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#FF3333", 0.4), medcol="#FF3333", border="#FF3333", outpch = 19, outcex = .5)

  aa <- mmM[[as.character(x)]][33:54,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#3333FF", 0.4), medcol="#3333FF", border="#3333FF", outpch = 19, outcex = .5)


  aa <- mmM[[as.character(x)]][1:3,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#466D1D", 0.4), medcol="#607D3B", border="#607D3B", outpch = 19, outcex = .5)
}

lo <- loess(unexp ~ sequ[-length(sequ)])
points(sequ[-length(sequ)], predict(lo), t = "l", col = "orange", lwd=3)
axis(1, at = seq(20, 1100, 4*15), seq(20, 1100, 4*15))
axis(2, at = seq(0, 1,0.1), paste0(seq(0,100,10), "%"), las=2)
axis(2, at = seq(0, 1,0.05), rep("", length(seq(0, 1,0.05))), las=2)
#legend("topleft", c("Male (G)", "Environment", "Male (GxG)"), fill = add.alpha(c("#FF3333", "#AFE1AF", "#3333FF"), 0.4), border = c("#FF3333", "#AFE1AF", "#3333FF"))

dev.off()



