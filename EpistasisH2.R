setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3")
regions <- read.table("regions_4way_merged_Oct24_effects.txt", sep = "\t", header=TRUE, row.names = 1)
regions[,"Top"] <- gsub(",", "", regions[,"Top"])

markers <- paste0(regions[, "Chr"], "_",regions[, "Top"], sep="")
names(markers) <- rownames(regions)

rmarkers <- gtsp[,unlist(lapply(markers, grep, colnames(gtsp)))]

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]),
                    GenoID = as.character(pull.pheno(mcross)[, "GenoID"]),
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

sequ <-seq(35, 1100, 15)
reps <- 1:50

all <- c("1_3010272", "1_24042124", "1_120474787", "2_89844287", "2_112712327", "2_148442635","3_83838529", "3_92135706", "4_52524395",
         "5_67573068", "6_107382038", "6_132762500", "9_29939029", "9_104091597", "9_124056586", "10_72780332", "11_5628810", "11_82178599",
         "12_112855820", "13_89689878", "14_101437457", "15_74248242", "17_32883804", "18_60822951")

names(all) <- c("Vita1a","Vita1b","Vita1c","Vita2a","Vita2b","Vita2c","Vita3a","Vita3b","Vita4a","Vita5a","Vita6a","Vita6b",
                "Vita9a","Vita9b","Vita9c","Vita10a","Vita11a","Vita11b","Vita12a","Vita13a","Vita14a","Vita15a",
                "Vita17a","Vita18a")
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3")
lodM.f <- read.table(paste0("vita_interactions_2way_females.txt"), sep = "\t")
lodM.m <- read.table(paste0("vita_interactions_2way_males.txt"), sep = "\t")

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
    Vita2a <- rmarkers[idx, 13:16]
    Vita2b <- rmarkers[idx, 17:20]
    Vita2c <- rmarkers[idx, 21:24]
    Vita3a <- rmarkers[idx, 25:28]
    Vita3b <- rmarkers[idx, 29:32]
    Vita4a <- rmarkers[idx, 33:36]
    Vita5a <- rmarkers[idx, 37:40]
    Vita6a <- rmarkers[idx, 41:44]
    Vita6b <- rmarkers[idx, 45:48]
    Vita9a <- rmarkers[idx, 49:52]
    Vita9b <- rmarkers[idx, 53:56]
    Vita9c <- rmarkers[idx, 57:60]
    Vita10a <- rmarkers[idx, 61:64]
    Vita11a <- rmarkers[idx, 65:68]
    Vita11b <- rmarkers[idx, 69:72]
    Vita12a <- rmarkers[idx, 73:76]
    Vita13a <- rmarkers[idx, 77:80]
    Vita14a <- rmarkers[idx, 81:84]
    Vita15a <- rmarkers[idx, 85:88]
    Vita17a <- rmarkers[idx, 89:92]
    Vita18a <- rmarkers[idx, 93:96]
    VitaXa <- rmarkers[idx, 97:100]
    VitaXb <- rmarkers[idx, 101:104]
    ml <- lm(longevity ~ site + cohort + treatment + 
             Vita1a + Vita1b + Vita1c + Vita2a + Vita2b + Vita2c + Vita3a + Vita3b + Vita4a + Vita5a + 
             Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita12a + Vita13a + 
             Vita14a + Vita15a + Vita17a + Vita18a + VitaXa + VitaXb + 
             Vita1a:Vita6a + Vita1a:Vita6b + Vita1a:Vita9a + Vita1a:Vita11a + Vita1b:Vita6b + Vita1b:Vita9a +
             Vita1b:Vita9b + Vita1c:Vita2b + Vita2b:Vita6a + Vita2b:Vita13a + Vita3a:Vita5a + Vita3b:Vita5a + 
             Vita3b:Vita12a + Vita4a:Vita13a + Vita5a:Vita17a + Vita6b:Vita14a, data = adata)
    ma <- anova(ml)
    # Mean square methods (most common one used)
    N <- (length(idx)/4)
    sQEE <- (ma[, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N            # Partial Total variance
    sQG <- (ma[4:45, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N         # Partial Genetic variance
    sQE <- (ma["Residuals", "Mean Sq"] + sum(ma[1:3, "Mean Sq"])) / N     # Fixed effects & residual variance
    sQTT <- (ma["Residuals", "Mean Sq"] + sum(ma[-46, "Mean Sq"])) / N    # Fixed effects & residual variance
    sQT <- sQG + sQE # Total variance
    vE <- sQEE / sQTT # Contributions
    Hs <- cbind(Hs, vE)
  }
  cat("Done ", x,"\n")
  colnames(Hs) <- reps
  rownames(Hs) <- c("site", "cohort", "treatment", c(names(markers), paste0("I", 1:16)), "unexplained")

  mmF[[as.character(x)]] <- Hs
}


mG <- c()
mE <- c()
for(x in sequ){
  aa <- round(apply(mmF[[as.character(x)]][4:45,],1, mean),3)
  mG <- rbind(mG, aa)
  aa <- round(apply(mmF[[as.character(x)]][1:3,],1, mean),3)
  mE <- rbind(mE, aa)
}
mG[mG < 0] <- 0
mE[mE < 0] <- 0
rownames(mG) <- sequ
colnames(mG) <- c(names(markers), paste0("I", 1:16))

rownames(mE) <- sequ
colnames(mE) <- c("Site", "Cohort", "Treatment")


add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny")
pdf("H2_G_GxG_Female.pdf", width = 26, height = 20)
op <- par(mar = c(4,10,2,2))
op <- par(cex = 1.5)

plot(c(10, 1110), y = c(0, 0.60), xlab = "Tage (Days)", ylab = "Heritability", 
     main = "Female Heritability - 26 Vita regions + 16 GxG interactions", t = "n", yaxt="n", xaxt="n", yaxs = "i", xaxs = "i")
unexp <- c()
for(x in 1:length(sequ)){
  unexp <- c(unexp, 1 - (sum(mG[x,]) + sum(mE[x,])))
}

for(x in sequ){
  aa <- mmF[[as.character(x)]][4:45,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#FF3333", 0.4), medcol="#FF3333", border="#FF3333", outpch = 19, outcex = .5)

  aa <- mmF[[as.character(x)]][30:45,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#3333FF", 0.4), medcol="#3333FF", border="#3333FF", outpch = 19, outcex = .5)


  aa <- mmF[[as.character(x)]][1:3,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#466D1D", 0.4), medcol="#607D3B", border="#607D3B", outpch = 19, outcex = .5)
}

lo <- loess(unexp ~ sequ)
points(sequ, predict(lo), t = "l")
axis(1, at = seq(20, 1100, 4*15), seq(20, 1100, 4*15))
axis(2, at = seq(0, 1,0.1), paste0(seq(0,100,10), "%"), las=2)
legend("topleft", c("Female (G + GxG)", "Environment", "Female (GxG)"), fill = add.alpha(c("#FF3333", "#AFE1AF", "#3333FF"), 0.4), border = c("#FF3333", "#AFE1AF", "#3333FF"))
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
    Vita2a <- rmarkers[idx, 13:16]
    Vita2b <- rmarkers[idx, 17:20]
    Vita2c <- rmarkers[idx, 21:24]
    Vita3a <- rmarkers[idx, 25:28]
    Vita3b <- rmarkers[idx, 29:32]
    Vita4a <- rmarkers[idx, 33:36]
    Vita5a <- rmarkers[idx, 37:40]
    Vita6a <- rmarkers[idx, 41:44]
    Vita6b <- rmarkers[idx, 45:48]
    Vita9a <- rmarkers[idx, 49:52]
    Vita9b <- rmarkers[idx, 53:56]
    Vita9c <- rmarkers[idx, 57:60]
    Vita10a <- rmarkers[idx, 61:64]
    Vita11a <- rmarkers[idx, 65:68]
    Vita11b <- rmarkers[idx, 69:72]
    Vita12a <- rmarkers[idx, 73:76]
    Vita13a <- rmarkers[idx, 77:80]
    Vita14a <- rmarkers[idx, 81:84]
    Vita15a <- rmarkers[idx, 85:88]
    Vita17a <- rmarkers[idx, 89:92]
    Vita18a <- rmarkers[idx, 93:96]
    VitaXa <- rmarkers[idx, 97:100]
    VitaXb <- rmarkers[idx, 101:104]
    ml <- lm(longevity ~ site + cohort + treatment + 
             Vita1a + Vita1b + Vita1c + Vita2a + Vita2b + Vita2c + Vita3a + Vita3b + Vita4a + Vita5a + 
             Vita6a + Vita6b + Vita9a + Vita9b + Vita9c + Vita10a + Vita11a + Vita11b + Vita12a + Vita13a + 
             Vita14a + Vita15a + Vita17a + Vita18a + VitaXa + VitaXb +
             Vita1c:Vita3a + Vita1c:Vita3b + Vita2a:Vita13a + Vita2b:Vita3a + Vita2b:Vita14a + Vita2c:Vita3b +
             Vita2c:Vita9b + Vita6b:Vita13a + Vita9b:Vita17a + Vita9c:Vita17a + Vita11b:Vita17a + Vita14a:Vita18a, data = adata);
    ma <- anova(ml)
    # Mean square methods (most common one used)
    N <- (length(idx)/4)
    sQEE <- (ma[, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N # Partial Total variance
    sQG <- (ma[4:41, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N # Partial Genetic variance
    sQE <- (ma["Residuals", "Mean Sq"] + sum(ma[1:3, "Mean Sq"])) / N # Fixed effects & residual variance
    sQTT <- (ma["Residuals", "Mean Sq"] + sum(ma[-42, "Mean Sq"])) / N # Fixed effects & residual variance
    sQT <- sQG + sQE # Total variance
    vE <- sQEE / sQTT # Contributions
    Hs <- cbind(Hs, vE)
  }
  colnames(Hs) <- reps
  rownames(Hs) <- c("site", "cohort", "treatment", c(names(markers), paste0("I", 1:12)), "unexplained")
  cat("Done ", x,"\n")
  mmM[[as.character(x)]] <- Hs
}



mG <- c()
mE <- c()
for(x in sequ){
  aa <- round(apply(mmM[[as.character(x)]][4:41,],1, mean),3)
  mG <- rbind(mG, aa)
  aa <- round(apply(mmM[[as.character(x)]][1:3,],1, mean),3)
  mE <- rbind(mE, aa)
}
mG[mG < 0] <- 0
mE[mE < 0] <- 0
rownames(mG) <- sequ
colnames(mG) <- c(names(markers), paste0("I", 1:12))

rownames(mE) <- sequ
colnames(mE) <- c("Site", "Cohort", "Treatment")



add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny")
pdf("H2_G_GxG_Male.pdf", width = 26, height = 20)
op <- par(mar = c(4,10,2,2))
op <- par(cex = 1.5)

plot(c(10, 1110), y = c(0, 0.60), xlab = "Tage (Days)", ylab = "Heritability", 
     main = "Male Heritability - 26 Vita regions + 12 GxG interactions", t = "n", yaxt="n", xaxt="n", yaxs = "i", xaxs = "i")
unexp <- c()
for(x in 1:length(sequ)){
  unexp <- c(unexp, 1 - (sum(mG[x,]) + sum(mE[x,])))
}

for(x in sequ){
  aa <- mmM[[as.character(x)]][4:41,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#FF3333", 0.4), medcol="#FF3333", border="#FF3333", outpch = 19, outcex = .5)

  aa <- mmM[[as.character(x)]][30:41,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#3333FF", 0.4), medcol="#3333FF", border="#3333FF", outpch = 19, outcex = .5)


  aa <- mmM[[as.character(x)]][1:3,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#466D1D", 0.4), medcol="#607D3B", border="#607D3B", outpch = 19, outcex = .5)
}

lo <- loess(unexp ~ sequ)
points(sequ, predict(lo), t = "l")
axis(1, at = seq(20, 1100, 4*15), seq(20, 1100, 4*15))
axis(2, at = seq(0, 1,0.1), paste0(seq(0,100,10), "%"), las=2)
legend("topleft", c("Male (G + GxG)", "Environment", "Male (GxG)"), fill = add.alpha(c("#FF3333", "#AFE1AF", "#3333FF"), 0.4), border = c("#FF3333", "#AFE1AF", "#3333FF"))

dev.off()
