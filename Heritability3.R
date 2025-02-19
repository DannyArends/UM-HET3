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
             Vita14a + Vita15a + Vita17a + Vita18a + VitaXa + VitaXb, data = adata)
    ma <- anova(ml)
    # Mean square methods (most common one used)
    N <- (length(idx)/4)
    sQEE <- (ma[, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N            # Partial Total variance
    sQG <- (ma[4:29, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N         # Partial Genetic variance
    sQE <- (ma["Residuals", "Mean Sq"] + sum(ma[1:3, "Mean Sq"])) / N     # Fixed effects & residual variance
    sQTT <- (ma["Residuals", "Mean Sq"] + sum(ma[-30, "Mean Sq"])) / N    # Fixed effects & residual variance
    sQT <- sQG + sQE # Total variance
    vE <- sQEE / sQTT # Contributions
    Hs <- cbind(Hs, vE)
  }
  cat("Done ", x,"\n")
  colnames(Hs) <- reps
  rownames(Hs) <- c("site", "cohort", "treatment", names(markers), "unexplained")

  mmF[[as.character(x)]] <- Hs
}


mG <- c()
mE <- c()
for(x in sequ){
  aa <- round(apply(mmF[[as.character(x)]][4:29,],1, mean),3)
  mG <- rbind(mG, aa)
  aa <- round(apply(mmF[[as.character(x)]][1:3,],1, mean),3)
  mE <- rbind(mE, aa)
}
mG[mG < 0] <- 0
mE[mE < 0] <- 0
rownames(mG) <- sequ
colnames(mG) <- c("Vita1a","Vita1b","Vita1c","Vita2a","Vita2b","Vita2c","Vita3a","Vita3b","Vita4a","Vita5a","Vita6a","Vita6b",
                "Vita9a","Vita9b","Vita9c","Vita10a","Vita11a","Vita11b","Vita12a","Vita13a","Vita14a","Vita15a",
                "Vita17a","Vita18a","VitaXa","VitaXb")

rownames(mE) <- sequ
colnames(mE) <- c("Site", "Cohort", "Treatment")

add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny")
pdf("H2_Female.pdf", width = 20, height = 20)
nf <- layout( matrix(c(1,2,2), ncol=1) )
op <- par(mar = c(4,10,2,2))
op <- par(cex = 1.5)
plot(c(10, 1110), y = c(0, 0.60), xlab = "Tage (Days)", ylab = "Heritability", 
     main = "Female Heritability - 26 Vita regions", t = "n", yaxt="n", xaxt="n", yaxs = "i", xaxs = "i")
for(x in sequ){
  aa <- mmF[[as.character(x)]][4:29,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#FF3333", 0.4), medcol="#FF3333", border="#FF3333", outpch = 19, outcex = .5)

  aa <- mmF[[as.character(x)]][1:3,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#466D1D", 0.4), medcol="#607D3B", border="#607D3B", outpch = 19, outcex = .5)
}
axis(1, at = seq(20, 1100, 4*15), seq(20, 1100, 4*15))
axis(2, at = seq(0, 1,0.1), paste0(seq(0,100,10), "%"), las=2)
legend("topleft", c("Female", "Environment"), fill = add.alpha(c("#FF3333", "#AFE1AF"), 0.4), border = c("#FF3333", "#AFE1AF"))


combined <- cbind(mE, round(100 * mG)[,ncol(mG):1])
ms <- c(c(0.0, 0.01, 0.02 ,0.04, 0.08, 0.16, 0.32, 0.64), c(1,2,4,8,16,32))
colz <- c(colorRampPalette(c("#FFFFFF","#7B57A4"))(7), colorRampPalette(c("#FFFFFF","#FF3333"))(6))
image(1:nrow(mG), 1:(ncol(mE) + ncol(mG)), combined, xaxt="n", 
      yaxt = "n", xlab = "Tage (Days)", ylab = "", main = "H2 Mean square methods (Females)",
      breaks = ms, col = colz)

axis(2, at = 1:ncol(combined), colnames(combined), las=2)
box()
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
             Vita14a + Vita15a + Vita17a + Vita18a + VitaXa + VitaXb, data = adata)
    ma <- anova(ml)
    # Variance explained: Nu - partial method
    #vE <- ma[, "Sum Sq"] / ma["Residuals", "Sum Sq"]


    # Variance explained: Omega partial method
    #vE <- (ma[, "Sum Sq"] - (ma[, "Df"] * ma["Residuals", "Mean Sq"])) / (ma[, "Sum Sq"] + ((length(idx) - ma[, "Df"]) * ma["Residuals", "Mean Sq"]))

    # Mean square methods (most common one used)
    N <- (length(idx)/4)
    sQEE <- (ma[, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N # Partial Total variance
    sQG <- (ma[4:29, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / N # Partial Genetic variance
    sQE <- (ma["Residuals", "Mean Sq"] + sum(ma[1:3, "Mean Sq"])) / N # Fixed effects & residual variance
    sQTT <- (ma["Residuals", "Mean Sq"] + sum(ma[-30, "Mean Sq"])) / N # Fixed effects & residual variance
    sQT <- sQG + sQE # Total variance
    vE <- sQEE / sQTT # Contributions
    Hs <- cbind(Hs, vE)
  }
  colnames(Hs) <- reps
  rownames(Hs) <- c("site", "cohort", "treatment", names(markers), "unexplained")
  cat("Done ", x,"\n")
  mmM[[as.character(x)]] <- Hs
}



mG <- c()
mE <- c()
for(x in sequ){
  aa <- round(apply(mmM[[as.character(x)]][4:29,],1, mean),3)
  mG <- rbind(mG, aa)
  aa <- round(apply(mmM[[as.character(x)]][1:3,],1, mean),3)
  mE <- rbind(mE, aa)
}
mG[mG < 0] <- 0
mE[mE < 0] <- 0
rownames(mG) <- sequ
colnames(mG) <- c("Vita1a","Vita1b","Vita1c","Vita2a","Vita2b","Vita2c","Vita3a","Vita3b","Vita4a","Vita5a","Vita6a","Vita6b",
                "Vita9a","Vita9b","Vita9c","Vita10a","Vita11a","Vita11b","Vita12a","Vita13a","Vita14a","Vita15a",
                "Vita17a","Vita18a","VitaXa","VitaXb")

rownames(mE) <- sequ
colnames(mE) <- c("Site", "Cohort", "Treatment")

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")
pdf("H2_Male.pdf", width = 20, height = 20)

nf <- layout( matrix(c(1,2,2), ncol=1) )
op <- par(mar = c(4,10,2,2))
op <- par(cex = 1.5)
plot(c(10, 1110), y = c(0, 0.70), xlab = "Tage (Days)", ylab = "Broad sense haplotype heritability", 
     main = "Male Heritability - 26 Vita regions", t = "n", yaxt="n", xaxt="n", yaxs = "i", xaxs = "i")
for(x in sequ){
  aa <- mmM[[as.character(x)]][4:29,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#00AEEF", 0.4), medcol="#00AEEF", border="#00AEEF", outpch = 19, outcex = .5)

  aa <- mmM[[as.character(x)]][1:3,]
  aa[aa < 0] <- 0
  boxplot(as.numeric(apply(aa,2,sum)), boxwex=20, at = x, add = TRUE, yaxt="n", 
          col = add.alpha("#466D1D", 0.4), medcol="#607D3B", border="#607D3B", outpch = 19, outcex = .5)
}
axis(1, at = seq(20, 1100, 4*15), seq(20, 1100, 4*15))
axis(2, at = seq(0, 1,0.1), paste0(seq(0,100,10), "%"), las=2)
legend("bottomleft", c("Male", "Environment"), fill = add.alpha(c("#00AEEF", "#AFE1AF"), 0.4), border = c("#00AEEF", "#AFE1AF"))

combined <- cbind(mE, round(100 * mG)[,ncol(mG):1])
ms <- c(c(0.0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64), c(1,2,4,8,16,32))
colz <- c(colorRampPalette(c("#FFFFFF","#486E33"))(7), colorRampPalette(c("#FFFFFF","#00AEEF"))(6))
image(1:nrow(mG), 1:(ncol(mE) + ncol(mG)), combined, xaxt="n", 
      yaxt = "n", xlab = "Tage (Days)", ylab = "", main = "H2 Mean square methods (Male)",
      breaks = ms, col = colz)

axis(2, at = 1:ncol(combined), colnames(combined), las=2)
box()
dev.off()


#### OLD


mR <- c()
for(x in sequ){
  x <- round(apply(mmM[[as.character(x)]][4:29,],1, mean),3)
  mR <- rbind(mR, x)
}
mR[mR < 0] <- 0
rownames(mR) <- sequ
colnames(mR) <- c("Vita1a","Vita1b","Vita1c","Vita2a","Vita2b","Vita2c","Vita3a","Vita3b","Vita4a","Vita5a","Vita6a","Vita6b",
                "Vita9a","Vita9b","Vita9c","Vita10a","Vita11a","Vita11b","Vita12a","Vita13a","Vita14a","Vita15a",
                "Vita17a","Vita18a","VitaXa","VitaXb")

#TODO from 42 days
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_Tables_Files/07_Heritability_Figures")
pdf("H2_Male.pdf", width = 20, height = 20)
  op <- par(mar = c(4,10,2,2))
  op <- par(cex = 1.5)
  image(1:nrow(mR), 1:ncol(mR), mR, xaxt="n", yaxt = "n", xlab = "Cut", ylab = "", main = "H2 Mean square methods (Males)")
  axis(1, at = (1:nrow(mR))[c(TRUE, FALSE)], sequ[c(TRUE, FALSE)])
  axis(2, at = 1:ncol(mR), colnames(mR), las=2)
  box()
dev.off()



#TODO from 42 days
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_Tables_Files/07_Heritability_Figures")
pdf(paste0("Extra_Heritability_Female_Env.pdf"), width = 14, height = 10)
  plot(c(365, 1100), y = c(0, 0.75), xlab = "Cutoff age", ylab = "Broad sense haplotype heritability", 
       main = "Heritability - 26 Vita regions", t = "n", yaxt="n", xaxt="n", yaxs = "i")
  for(x in sequ){
    vita <- apply(mmF[[as.character(x)]],2, function(y){sum(y[4:29])})
    boxplot(as.numeric(vita), boxwex=20, at = x-2.5, add = TRUE, yaxt="n", col = "hotpink")

    env <- apply(mmF[[as.character(x)]],2, function(y){sum(y[1:3])})
    boxplot(env, boxwex=10, at = x+2.5, add = TRUE, yaxt="n", col = 3)

  #  vita <- apply(mmM[[as.character(x)]],2, function(y){sum(y[4:29])})
  #  boxplot(as.numeric(vita), boxwex=20, at = x-2.5, add = TRUE, yaxt="n", col = "blue")

  #  env <- apply(mmM[[as.character(x)]],2, function(y){sum(y[1:3])})
  #  boxplot(env, boxwex=10, at = x+2.5, add = TRUE, yaxt="n", col = 3)

    cat(mean(vita), "\n")
  }
  axis(1, at = seq(365, 1100, 4*15), seq(365, 1100, 4*15))
  axis(2, at = seq(0, 1,0.1), paste0(seq(0,100,10), "%"), las=2)
dev.off()


