library(svglite)

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
                    adjLongevity = NA, 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw6 = as.numeric(pull.pheno(mcross)[, "bw6"]),
                    adjBw6 = NA,
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
#idx <- which(cdata[, "longevity"] >= 365 & cdata[, "sex"] == 0) # females
#idx <- which(cdata[, "longevity"] >= 365 & cdata[, "sex"] == 1) # males
idx <- which(cdata[, "longevity"] >= 365)
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
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
#pdf(paste0("Figure_4_Vita9C_Combined.pdf"), width = 8, height = 8)

plot(c(365, 1250), c(0, 100), t = "n", ylab = "% survival", xlab = "days", yaxt="n", main = "KM Curve (Vita9C - ALL)")
for(x in 365:1456) {
  n.C3 <- length(which(cdata[C3, "adjLongevity"] >= x))
  n.D2 <- length(which(cdata[D2, "adjLongevity"]  >= x))

  points(x, (n.C3/length(C3))* 100, pch=19, cex=0.2, col = "red")
  points(x, (n.D2/length(D2))* 100, pch=19, cex=0.2, col = "green")
}
axis(2, at = c(0,25,50,75,100), c(0,25,50,75,100), las=2)
legend("topright", c("C3", "D2"), col = c("red", "green"), pch=18)
#dev.off()

### 4 Haplotypes

CH <- which(cdata[, "gts"] == "AC")
CD <- which(cdata[, "gts"] == "AD")
BH <- which(cdata[, "gts"] == "BC")
BD <- which(cdata[, "gts"] == "BD")

col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
col.main <- adjustcolor( col.main, alpha.f = 0.6)

plot(c(365, 1250), c(0, 100), t = "n", ylab = "% survival", xlab = "days", yaxt="n", main = "KM Curve")
for(x in 365:1456) {
  n.CH <- length(which(cdata[CH, "adjLongevity"] >= x))
  n.CD <- length(which(cdata[CD, "adjLongevity"]  >= x))
  n.BH <- length(which(cdata[BH, "adjLongevity"]  >= x))
  n.BD <- length(which(cdata[BD, "adjLongevity"]  >= x))

  points(x, (n.CH/length(CH))* 100, pch=19, cex=0.2, col = col.main[1])
  points(x, (n.CD/length(CD))* 100, pch=19, cex=0.2, col = col.main[2])
  points(x, (n.BH/length(BH))* 100, pch=19, cex=0.2, col = col.main[3])
  points(x, (n.BD/length(BD))* 100, pch=19, cex=0.2, col = col.main[4])

}
axis(2, at = c(0,25,50,75,100), c(0,25,50,75,100), las=2)
legend("topright", c("CH", "CD", "BH", "BD"), col = col.main, pch=18)

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
  means <- c(mean(pheAdj[which(gts == "AC" | gts== "BC")]),mean(pheAdj[which(gts == "AD" | gts == "DB")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "AC"| gts== "BC")]),std(pheAdj[which(gts == "AD"| gts == "DB")]))
  return(list(
          c(round(OAmean,0), round(means,0)),
          c(round(stderrs,2))
        ))
}

pos <- "9_104091597"
gtsp <- pull.genoprob(mcross)

remaining <- c()
errors <- c()
msequence <- seq(20, 1100, 15)
for(d in msequence){
  combined <- getEffect(mcross, gtsp, marker = pos, timepoint = d)
  female <- getEffect(mcross, gtsp, marker = pos, timepoint = d, sex = 0, model = "longevity ~ site + cohort + treatment")
  male <- getEffect(mcross, gtsp, marker = pos, timepoint = d, sex = 1, model = "longevity ~ site + cohort + treatment")
  remaining <- rbind(remaining, c(combined[[1]], female[[1]], male[[1]]))
  errors <- rbind(errors, c(combined[[2]], female[[2]], male[[2]]))
}

col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
col.alpha <- add.alpha(col.main, 0.1)

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/11_FiguresRedone")
pdf(paste0("Vita9B_C3_D2_alleleEffects.pdf"), width = 36, height = 12)
  op <- par(mfrow = c(1,3))
  op <- par(cex = 2)
  par(mar = c(5, 5, 4, 3))

  plot(c(20, 1100), c(-25, 25), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
         ylab = "Actuarial Effect Size [d]", main = paste0("Vita9b Combined"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
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
  legend("top", c("C3", "D2"), col = col.main[1:2], lwd=4, bg = "white", ncol=2, bty = "n")
  box()

  # Females
  plot(c(20, 1100), c(-25, 25), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
        ylab = "", main = paste0("Vita9b Females"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2, cex.axis = 1.4)

  points(msequence, remaining[,5], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,5] + errors[,3], rev(remaining[,5] - errors[,3])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,6], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,6] + errors[,4], rev(remaining[,6] - errors[,4])), col = col.alpha[2], border = NA)
  legend("top", c("C3", "D2"), col = col.main[1:2], lwd=4, bg = "white", ncol=2, bty = "n")
  box()

  # Males
  plot(c(20, 1100), c(-25, 25), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
         ylab = "", main = paste0("Vita9b Males"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2,cex.axis = 1.4)

  points(msequence, remaining[,8], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,8] + errors[,5], rev(remaining[,8] - errors[,5])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,9], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,9] + errors[,6], rev(remaining[,9] - errors[,6])), col = col.alpha[2], border = NA)

  legend("top", c("C3", "D2"), col = col.main[1:2], lwd=4, bg = "white", ncol=2, bty = "n")
  box()
#dev.off()
