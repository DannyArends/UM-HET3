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

se <- sqrt(cumsum(onTP["n.event"] / (onTP["n.risk"] *(onTP["n.risk"] - onTP["n.event"]))))

tp <- sort(unique(cdata[C3, "adjLongevity"]))
res.C3 <- c()
n.evt <- 0
for(x in tp) {
  n.C3 <- length(which(cdata[C3, "adjLongevity"] == x))
  onTP <- c(time = x, n.risk =  length(C3) - n.evt, n.event = n.C3, survival = (length(C3) - (n.evt + n.C3)) / length(C3))
  res.C3 <- rbind(res.C3, onTP)
  n.evt <- n.evt + n.C3
}
std <- sqrt(cumsum(res.C3[,"n.event"] / (res.C3[,"n.risk"] *(res.C3[,"n.risk"] - res.C3[,"n.event"]))))
se <- std * res.C3[, "survival"]
cbind(res.C3, std.err = se)

chaz <- -log(res.C3[,"survival"])
zstat <- -qnorm((1 - 0.90)/2)
lower <- exp(-(chaz + zstat*std))
upper <- exp(-(chaz - zstat*std))
plot(tp, upper, t = "l", col = "red")
points(tp, lower, t = "l", col = "green")

sqrt(cumsum(res.C3$n.event / (res.C3$n.risk * (res.C3$n.risk - res.C3$n.event))))

#plot(c(365, 1250), c(0, 100), t = "n", ylab = "% survival", xlab = "days", yaxt="n", main = "KM Curve (Vita9C - ALL)")
for(x in 365:1456) {
  n.C3 <- length(which(cdata[C3, "adjLongevity"] >= x))
  n.D2 <- length(which(cdata[D2, "adjLongevity"]  >= x))

  points(x, (n.C3/length(C3))* 1, pch=19, cex=0.2, col = "purple")
 # points(x, (n.D2/length(D2))* 100, pch=19, cex=0.2, col = "green")
}
axis(2, at = c(0,25,50,75,100), c(0,25,50,75,100), las=2)
legend("topright", c("C3", "D2"), col = c("red", "green"), pch=18)
