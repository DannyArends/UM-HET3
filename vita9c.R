setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)

vita9c <- "9_104091597"

mm <- pull.geno(fill.geno(mcross))
phe <- pull.pheno(mcross)[, "longevity"]

sex <- pull.pheno(mcross)[, "sex"]
site <- as.factor(pull.pheno(mcross)[, "site"])
cohort <- as.factor(pull.pheno(mcross)[, "cohort"])
treatment <- as.factor(pull.pheno(mcross)[, "treatment"])

marker <- mm[,vita9c]

Wf <- which(sex == 0)
Wm <- which(sex == 1)

curves <- c()
msequence <- seq(1, max(phe), 1)
for(d in msequence){
  all = round(100 * length(which(phe > d)) / length(phe),0)
  females = round(100 * length(which(phe[Wf] > d)) / length(phe[Wf]),1)
  males = round(100 * length(which(phe[Wm] > d)) / length(phe[Wm]),1)

  split.f <- c()
  split.m <- c()
  for(x in 1:4){
    If <- which(sex == 0 & marker == x)
    Im <- which(sex == 1 & marker == x)
    split.f <- c(split.f, round(100 * length(which(phe[If] > d)) / length(phe[If]),1))
    split.m <- c(split.m, round(100 * length(which(phe[Im] > d)) / length(phe[Im]),1))
  }

  curves <- rbind(curves, c(all, females, split.f, males, split.m))
}

curves[curves < 1] <- NA
curves <- log10(curves)

plot(c(500, 1300), c(-2, 2), t = 'n', xlab = "days", ylab = "log10(% surviving)", main = paste0("Survival curve at vita9c"))
rect(0,0,2000, 3, col = rgb(1,0.8,0.8,0.2))
rect(0,0,2000, -3, col = rgb(0,0.0,0.8,0.2))

# Females
points(msequence, curves[,2], t = 'l', col = "black", lwd=2)
points(msequence, curves[,3], t = 'l', col = "red", lwd=2)
points(msequence, curves[,4], t = 'l', col = "green", lwd=2)
points(msequence, curves[,5], t = 'l', col = "blue", lwd=2)
points(msequence, curves[,6], t = 'l', col = "orange", lwd=2)

# Males
points(msequence, -curves[,7], t = 'l', col = "black", lwd=2)
points(msequence, -curves[,8], t = 'l', col = "red", lwd=2)
points(msequence, -curves[,9], t = 'l', col = "green", lwd=2)
points(msequence, -curves[,10], t = 'l', col = "blue", lwd=2)
points(msequence, -curves[,11], t = 'l', col = "orange", lwd=2)

legend("topright", c("Avg", "C||H", "C||D", "B||H", "B||D"), col = c("black", "red", "green", "blue", "orange"), lwd=2)



remaining <- c()
msequence <- seq(365, 1100, 15)
for(d in msequence){
  all = mean(phe[which(phe > d)])

  females = mean(phe[Wf[which(phe[Wf] > d)]])
  males = mean(phe[Wm[which(phe[Wm] > d)]])

  rem.f <- c()
  rem.m <- c()
  for(x in 1:4){
    If <- which(sex == 0 & marker == x)
    Im <- which(sex == 1 & marker == x)
    rem.f <- c(rem.f, mean(phe[If[which(phe[If] > d)]]))
    rem.m <- c(rem.m, mean(phe[Im[which(phe[Im] > d)]]))
  }
  remaining <- rbind(remaining, c(all, females, rem.f, males, rem.m))
}

remaining <- round(remaining,0)
remaining[,3] <- remaining[,3] - remaining[,2]
remaining[,4] <- remaining[,4] - remaining[,2]
remaining[,5] <- remaining[,5] - remaining[,2]
remaining[,6] <- remaining[,6] - remaining[,2]

remaining[,8] <- remaining[,8] - remaining[,7]
remaining[,9] <- remaining[,9] - remaining[,7]
remaining[,10] <- remaining[,10] - remaining[,7]
remaining[,11] <- remaining[,11] - remaining[,7]

op <- par(mfrow = c(1,2))
plot(c(365, 1100), c(-20, 20), t = 'n', xlab = "days", ylab = "Allele effect (days)", main = paste0("Allele effect at vita9c (females)"))
# Females
points(msequence, remaining[,3], t = 'l', col = "red", lwd=2)
points(msequence, remaining[,4], t = 'l', col = "green", lwd=2)
points(msequence, remaining[,5], t = 'l', col = "blue", lwd=2)
points(msequence, remaining[,6], t = 'l', col = "orange", lwd=2)
legend("topleft", c("C||H", "C||D", "B||H", "B||D"), col = c("red", "green", "blue", "orange"), lwd=2)

# Males
plot(c(365, 1100), c(-20, 20), t = 'n', xlab = "days", ylab = "Allele effect (days)", main = paste0("Allele effect at vita9c (males)"))
points(msequence, remaining[,8], t = 'l', col = "red", lwd=2)
points(msequence, remaining[,9], t = 'l', col = "green", lwd=2)
points(msequence, remaining[,10], t = 'l', col = "blue", lwd=2)
points(msequence, remaining[,11], t = 'l', col = "orange", lwd=2)

legend("topleft", c("C||H", "C||D", "B||H", "B||D"), col = c("red", "green", "blue", "orange"), lwd=2)


