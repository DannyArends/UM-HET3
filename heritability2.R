setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

regions <- read.table("regions_4way_sorted_final.txt", sep = "\t", header=TRUE)
regions[,"Top"] <- gsub(",", "", regions[,"Top"])

markers <- paste0(regions[, "Chr"], "_",regions[, "Top"], sep="")
names(markers) <- regions[,1]

rmarkers <- gtsp[,unlist(lapply(markers, grep, colnames(gtsp)))]

cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]),
                    GenoID = as.character(pull.pheno(mcross)[, "GenoID"]),
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

sequ <-seq(365, 1100, 15)
reps <- 1:10

mm <- vector("list", length(sequ))
names(mm) <- sequ
for(x in sequ){
  Hs <- c()
  for(n in reps){
    idx <- which(cdata[, "longevity"] >= x)
    idx <- sample(idx, length(idx) * .7)
    adata <- cdata[idx, ]
    Vita1a <- rmarkers[idx, 1:4]
    Vita1b <- rmarkers[idx, 5:8]
    Vita1c <- rmarkers[idx, 9:12]
    Vita2a <- rmarkers[idx, 13:16]
    Vita2b <- rmarkers[idx, 17:20]
    Vita3a <- rmarkers[idx, 21:24]
    Vita3b <- rmarkers[idx, 25:28]
    Vita4a <- rmarkers[idx, 29:32]
    Vita5a <- rmarkers[idx, 33:36]
    Vita6a <- rmarkers[idx, 37:40]
    Vita6b <- rmarkers[idx, 41:44]
    Vita9a <- rmarkers[idx, 45:48]
    Vita9c <- rmarkers[idx, 49:52]
    Vita9d <- rmarkers[idx, 53:56]
    Vita10a <- rmarkers[idx, 57:60]
    Vita11a <- rmarkers[idx, 61:64]
    Vita11b <- rmarkers[idx, 65:68]
    Vita12a <- rmarkers[idx, 69:72]
    Vita13a <- rmarkers[idx, 73:76]
    Vita14a <- rmarkers[idx, 77:80]
    Vita15a <- rmarkers[idx, 81:84]
    Vita17a <- rmarkers[idx, 85:88]
    Vita17b <- rmarkers[idx, 89:92]
    Vita18a <- rmarkers[idx, 93:96]
    VitaXa <- rmarkers[idx, 97:100]
    VitaXb <- rmarkers[idx, 101:104]
    ml <- lm(longevity ~ sex + site + cohort + treatment + 
             Vita1a + Vita1b + Vita1c + Vita2a + Vita2b + Vita3a + Vita3b + Vita4a + Vita5a + Vita6a + Vita6b +
             Vita9a + Vita9c + Vita9d + Vita10a + Vita11a + Vita11b + Vita12a +Vita13a + Vita14a + Vita15a + 
             Vita17a + Vita17b + Vita18a + VitaXa + VitaXb, data = adata)
    ma <- anova(ml)
    # Variance explained: Nu - partial method
    vE <- ma[, "Sum Sq"] / ma["Residuals", "Sum Sq"]


    # Variance explained: Omega partial method
    #vE <- (ma[, "Sum Sq"] - (ma[, "Df"] * ma["Residuals", "Mean Sq"])) / (ma[, "Sum Sq"] + ((length(idx) - ma[, "Df"]) * ma["Residuals", "Mean Sq"]))

    # Mean square methods (most common one used)
    #sQG <- (ma[, "Mean Sq"] - ma["Residuals", "Mean Sq"]) / (length(idx)/4) # Partial Genetic variance
    #sQE <- (ma["Residuals", "Mean Sq"] + sum(ma[1:4, "Mean Sq"])) / (length(idx)/4) # Fixed effects & residual variance
    #sQT <- sQG + sQE # Total variance
    #vE <- sQG / sQT # Contributions
    Hs <- cbind(Hs, vE)
  }
  colnames(Hs) <- reps
  rownames(Hs) <- c("sex", "site", "cohort", "treatment", names(markers), "unexplained")

  mm[[as.character(x)]] <- Hs
}

plot(c(365, 1100), y = c(0, .4), xlab = "Cutoff age", ylab = "Broad sense haplotype heritability", main = "Combined 26 Vita regions", t = "n")
for(x in sequ){
  vita <- apply(mm[[as.character(x)]],2, function(y){sum(y[5:30])})
  #vita <- mm[[as.character(x)]]
  boxplot(as.numeric(vita), boxwex=10, at = x-2.5, add = TRUE, yaxt="n")

  env <- apply(mm[[as.character(x)]],2, function(y){sum(y[1:4])})
  boxplot(env, boxwex=10, at = x+2.5, add = TRUE, yaxt="n", col = "green")
  cat(mean(vita), "\n")
}

