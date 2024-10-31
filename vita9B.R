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
pdf(paste0("Figure_4_Vita9C_Combined.pdf"), width = 8, height = 8)

plot(c(365, 1250), c(0, 100), t = "n", ylab = "% survival", xlab = "days", yaxt="n", main = "KM Curve (Vita9C - ALL)")
for(x in 365:1456) {
  n.C3 <- length(which(cdata[C3, "adjLongevity"] >= x))
  n.D2 <- length(which(cdata[D2, "adjLongevity"]  >= x))

  points(x, (n.C3/length(C3))* 100, pch=19, cex=0.2, col = "red")
  points(x, (n.D2/length(D2))* 100, pch=19, cex=0.2, col = "green")
}
axis(2, at = c(0,25,50,75,100), c(0,25,50,75,100), las=2)
legend("topright", c("C3", "D2"), col = c("red", "green"), pch=18)
dev.off()

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

