#
# Vita9B.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# In Depth analysis and visualization of the effect of Vita9B, comparing back to a previously detected QTL 
# The visualization is seen in Figure 6a,b,c,d
#

library(svglite)

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny/Haplotypes_All")

gtsp <- pull.genoprob(mcross)
cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    adjLongevity = NA, 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw6 = as.numeric(pull.pheno(mcross)[, "bw6"]),
                    adjBw6 = NA,
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
#idx <- which(cdata[, "longevity"] >= 365 & cdata[, "sex"] == 0) # females
#idx <- which(cdata[, "longevity"] >= 365 & cdata[, "sex"] == 1) # males
idx <- which(cdata[, "longevity"] >= 20)
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),0)

mp <- gtsp[, grep("9_104091597", colnames(gtsp))]
gts <- unlist(lapply(lapply(lapply(apply(mp,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
  if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
}))

cdata <- cbind(cdata, gts)
cdata <- cdata[!is.na(cdata[, "gts"]),]
cdata <- cdata[which(cdata[, "adjLongevity"] >= 365),]

CH <- which(cdata[, "gts"] == "AC")
CD <- which(cdata[, "gts"] == "AD")
BH <- which(cdata[, "gts"] == "BC")
BD <- which(cdata[, "gts"] == "BD")

C3 <- c(CH, BH)
D2 <- c(CD, BD)
#setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresDanny")
#pdf(paste0("Figure_4_Vita9B_Combined.pdf"), width = 16, height = 8)
#op <- par(mfrow = c(1,2))
#plot(c(365, 1250), c(0, 100), t = "n", ylab = "% survival", xlab = "days", yaxt="n", main = "KM Curve (Vita9B - ALL)")
#for(x in 365:1456) {
#  n.C3 <- length(which(cdata[C3, "adjLongevity"] >= x))
#  n.D2 <- length(which(cdata[D2, "adjLongevity"]  >= x))

#  points(x, (n.C3/length(C3))* 100, pch=19, cex=0.2, col = "#00A654")
#  points(x, (n.D2/length(D2))* 100, pch=19, cex=0.2, col = "#004BAD")
#}
#axis(2, at = c(0,25,50,75,100), c(0,25,50,75,100), las=2)
#legend("topright", c("C3", "D2"), col = c("#00A654", "#004BAD"), pch=18)


### 4 Haplotypes

#CH <- which(cdata[, "gts"] == "AC")
#CD <- which(cdata[, "gts"] == "AD")
#BH <- which(cdata[, "gts"] == "BC")
#BD <- which(cdata[, "gts"] == "BD")

#col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
#col.main <- adjustcolor( col.main, alpha.f = 0.6)

#plot(c(365, 1250), c(0, 100), t = "n", ylab = "% survival", xlab = "days", yaxt="n", main = "KM Curve")
#for(x in 365:1456) {
#  n.CH <- length(which(cdata[CH, "adjLongevity"] >= x))
#  n.CD <- length(which(cdata[CD, "adjLongevity"]  >= x))
#  n.BH <- length(which(cdata[BH, "adjLongevity"]  >= x))
#  n.BD <- length(which(cdata[BD, "adjLongevity"]  >= x))

#  points(x, (n.CH/length(CH))* 100, pch=19, cex=0.2, col = col.main[1])
#  points(x, (n.CD/length(CD))* 100, pch=19, cex=0.2, col = col.main[2])
#  points(x, (n.BH/length(BH))* 100, pch=19, cex=0.2, col = col.main[3])
#  points(x, (n.BD/length(BD))* 100, pch=19, cex=0.2, col = col.main[4])

#}
#axis(2, at = c(0,25,50,75,100), c(0,25,50,75,100), las=2)
#legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=18)
#dev.off()

getEffect <- function(mcross, gtsprob, timepoint = 365, sex = "all", marker = "1_3010274", model = "longevity ~ sex + site + cohort + treatment"){
  mp <- gtsprob[, grep(marker, colnames(gtsprob))]
  gts <- unlist(lapply(lapply(lapply(apply(mp,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]),
                      geno = gts)
  idx <- which(cdata[, "longevity"] >= timepoint)
  cdata <- cdata[idx,]
  if(sex == 0 || sex == 1){ # males = 1, fem = 0
    cdata <- cdata[which(cdata[, "sex"] == sex),]
  }
  #cat("N=", nrow(cdata), "\n")

  gts <- cdata[,"geno"]
  mlm <- lm(as.formula(model), data = cdata)

  pheAdj <- rep(NA, nrow(cdata))
  names(pheAdj) <-  rownames(cdata)

  adj <- residuals(mlm) + mean(cdata[, "longevity"])
  pheAdj[names(adj)] <- adj
  OAmean <- mean(pheAdj[which(!is.na(gts))])
  means <- c(mean(pheAdj[which(gts == "AC" | gts== "BC")]),mean(pheAdj[which(gts == "AD" | gts == "BD")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "AC"| gts== "BC")]),std(pheAdj[which(gts == "AD"| gts == "BD")]))
  test <- t.test(pheAdj[which(gts == "AC" | gts== "BC")], pheAdj[which(gts == "AD" | gts == "BD")])
  return(list(
          c(round(OAmean,0), round(means,0)),
          c(round(stderrs,2)),
          test
        ))
}

all <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085")

names(all) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa")

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__bioRxiv_All_Key_Files/11_FiguresDanny/AlleleEffects")

for(pname in names(all)){
pos <- all[pname]#"9_104091597"
gtsp <- pull.genoprob(mcross)

remaining <- c()
errors <- c()
lods <- c()
msequence <- seq(20, 1100, 15)
for(d in msequence){
  combined <- getEffect(mcross, gtsp, marker = pos, timepoint = d)
  female <- getEffect(mcross, gtsp, marker = pos, timepoint = d, sex = 0, model = "longevity ~ site + cohort + treatment")
  male <- getEffect(mcross, gtsp, marker = pos, timepoint = d, sex = 1, model = "longevity ~ site + cohort + treatment")
  remaining <- rbind(remaining, c(combined[[1]], female[[1]], male[[1]]))
  errors <- rbind(errors, c(combined[[2]], female[[2]], male[[2]]))
  lods <- rbind(lods, c(combined[[3]]$p.value, female[[3]]$p.value, male[[3]]$p.value))
}

col.main <- c("#ff1493", "#deb887", "#B16BE6", "#F02D27")
add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
col.alpha <- add.alpha(col.main, 0.1)


pdf(paste0(pname,"_H_D_alleleEffects.pdf"), width = 36, height = 12)
  op <- par(mfrow = c(1,3))
  op <- par(cex = 2)
  par(mar = c(5, 5, 4, 3))

  plot(c(20, 1100), c(-25, 25), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
         ylab = "Actuarial Effect Size [d]", main = paste0(pname," Combined"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2,cex.axis = 1.4)
  #abline(v = seq(400, 1100, 50), lty=2, col="lightgray")
  #abline(h = seq(-50, 50, 10), lty=2, col="lightgray")

  # Combined
  points(msequence, remaining[,2], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,2] + errors[,1], rev(remaining[,2] - errors[,1])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,3], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,3] + errors[,2], rev(remaining[,3] - errors[,2])), col = col.alpha[2], border = NA)
  legend("top", c("H", "D"), col = col.main[1:2], lwd=4, bg = "white", ncol=2, bty = "n")

  abline(h = -25 + 3.95, lty = 1, col = "orange")
  points(seq(20, 1100, 15), -log10(lods[,1]) - 25, t = "l", lwd=2)
  axis(4, at = seq(-40, -32, 2), seq(0, 8, 2), las = 2,  cex.axis =1.1)
  box()

  # Females
  plot(c(20, 1100), c(-25, 25), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
        ylab = "", main = paste0(pname," Females"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2, cex.axis = 1.4)

  points(msequence, remaining[,5], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,5] + errors[,3], rev(remaining[,5] - errors[,3])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,6], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,6] + errors[,4], rev(remaining[,6] - errors[,4])), col = col.alpha[2], border = NA)
  legend("top", c("H", "D"), col = col.main[1:2], lwd=4, bg = "white", ncol=2, bty = "n")

  abline(h = -25 + 3.95, lty = 1, col = "orange")
  points(seq(20, 1100, 15), -log10(lods[,2]) - 25, t = "l", lwd=2)
  axis(4, at = seq(-40, -32, 2), seq(0, 8, 2), las = 2,  cex.axis =1.1)
  box()

  # Males
  plot(c(20, 1100), c(-25, 25), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
         ylab = "", main = paste0(pname," Males"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2,cex.axis = 1.4)

  points(msequence, remaining[,8], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,8] + errors[,5], rev(remaining[,8] - errors[,5])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,9], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,9] + errors[,6], rev(remaining[,9] - errors[,6])), col = col.alpha[2], border = NA)

  legend("top", c("H", "D"), col = col.main[1:2], lwd=4, bg = "white", ncol=2, bty = "n")

  abline(h = -25 + 3.95, lty = 1, col = "orange")
  points(seq(20, 1100, 15), -log10(lods[,3]) - 25, t = "l", lwd=2)
  axis(4, at = seq(-40, -32, 2), seq(0, 8, 2), las = 2,  cex.axis =1.1)
  box()
dev.off()
cat("Done", pname, "\n")
}

getEffect <- function(mcross, gtsprob, timepoint = 365, sex = "all", marker = "1_3010274", model = "longevity ~ sex + site + cohort + treatment"){
  mp <- gtsprob[, grep(marker, colnames(gtsprob))]
  gts <- unlist(lapply(lapply(lapply(apply(mp,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]),
                      geno = gts)
  idx <- which(cdata[, "longevity"] >= timepoint)
  cdata <- cdata[idx,]
  if(sex == 0 || sex == 1){ # males = 1, fem = 0
    cdata <- cdata[which(cdata[, "sex"] == sex),]
  }
  cat("N=", nrow(cdata), "\n")

  gts <- cdata[,"geno"]
  mlm <- lm(as.formula(model), data = cdata)

  pheAdj <- rep(NA, nrow(cdata))
  names(pheAdj) <-  rownames(cdata)

  adj <- residuals(mlm) + mean(cdata[, "longevity"])
  pheAdj[names(adj)] <- adj
  OAmean <- mean(pheAdj[which(!is.na(gts))])
  means <- c(mean(pheAdj[which(gts == "AC" | gts== "AD")]),mean(pheAdj[which(gts == "BC" | gts == "BD")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "AC"| gts== "AD")]),std(pheAdj[which(gts == "BC"| gts == "BD")]))
  test <- t.test(pheAdj[which(gts == "AC" | gts== "AD")], pheAdj[which(gts == "BC" | gts == "BD")])
  return(list(
          c(round(OAmean,0), round(means,0)),
          c(round(stderrs,2)),
          test
        ))
}

for(pname in names(all)[29]){
pos <- all[pname]#"9_104091597"
gtsp <- pull.genoprob(mcross)

remaining <- c()
errors <- c()
lods <- c()
msequence <- seq(20, 1100, 15)
for(d in msequence){
  combined <- getEffect(mcross, gtsp, marker = pos, timepoint = d)
  female <- getEffect(mcross, gtsp, marker = pos, timepoint = d, sex = 0, model = "longevity ~ site + cohort + treatment")
  male <- getEffect(mcross, gtsp, marker = pos, timepoint = d, sex = 1, model = "longevity ~ site + cohort + treatment")
  remaining <- rbind(remaining, c(combined[[1]], female[[1]], male[[1]]))
  errors <- rbind(errors, c(combined[[2]], female[[2]], male[[2]]))
  lods <- rbind(lods, c(combined[[3]]$p.value, female[[3]]$p.value, male[[3]]$p.value))
}

col.main <- c("#00cdcd", "#9400d3", "#B16BE6", "#F02D27")
add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
col.alpha <- add.alpha(col.main, 0.1)


pdf(paste0(pname,"_C_B_alleleEffects.pdf"), width = 36, height = 12)
  op <- par(mfrow = c(1,3))
  op <- par(cex = 2)
  par(mar = c(5, 5, 4, 3))

  plot(c(20, 1100), c(-25, 25), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
         ylab = "Actuarial Effect Size [d]", main = paste0(pname," Combined"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2,cex.axis = 1.4)
  #abline(v = seq(400, 1100, 50), lty=2, col="lightgray")
  #abline(h = seq(-50, 50, 10), lty=2, col="lightgray")

  # Combined
  points(msequence, remaining[,2], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,2] + errors[,1], rev(remaining[,2] - errors[,1])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,3], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,3] + errors[,2], rev(remaining[,3] - errors[,2])), col = col.alpha[2], border = NA)
  legend("top", c("C", "B"), col = col.main[1:2], lwd=4, bg = "white", ncol=2, bty = "n")

  abline(h = -25 + 3.95, lty = 1, col = "orange")
  points(seq(20, 1100, 15), -log10(lods[,1]) - 25, t = "l", lwd=2)
  axis(4, at = seq(-40, -32, 2), seq(0, 8, 2), las = 2,  cex.axis =1.1)
  box()

  # Females
  plot(c(20, 1100), c(-25, 25), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
        ylab = "", main = paste0(pname," Females"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2, cex.axis = 1.4)

  points(msequence, remaining[,5], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,5] + errors[,3], rev(remaining[,5] - errors[,3])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,6], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,6] + errors[,4], rev(remaining[,6] - errors[,4])), col = col.alpha[2], border = NA)
  legend("top", c("C", "B"), col = col.main[1:2], lwd=4, bg = "white", ncol=2, bty = "n")

  abline(h = -25 + 3.95, lty = 1, col = "orange")
  points(seq(20, 1100, 15), -log10(lods[,2]) - 25, t = "l", lwd=2)
  axis(4, at = seq(-40, -32, 2), seq(0, 8, 2), las = 2,  cex.axis =1.1)
  box()

  # Males
  plot(c(20, 1100), c(-25, 25), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
         ylab = "", main = paste0(pname," Males"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2,cex.axis = 1.4)

  points(msequence, remaining[,8], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,8] + errors[,5], rev(remaining[,8] - errors[,5])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,9], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,9] + errors[,6], rev(remaining[,9] - errors[,6])), col = col.alpha[2], border = NA)

  legend("top", c("C", "B"), col = col.main[1:2], lwd=4, bg = "white", ncol=2, bty = "n")

  abline(h = -25 + 3.95, lty = 1, col = "orange")
  points(seq(20, 1100, 15), -log10(lods[,3]) - 25, t = "l", lwd=2)
  axis(4, at = seq(-40, -32, 2), seq(0, 8, 2), las = 2,  cex.axis =1.1)
  box()
dev.off()
}

