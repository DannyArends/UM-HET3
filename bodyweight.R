setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

# Our Progressive Mapping Sequence
msequence <- seq(365, 1100, 15)
markers <- unique(unlist(lapply(strsplit(colnames(gtsp), ":"), "[",1)))

lods.cM <- c()
minAge <- c()
cdata <- data.frame(bw6 = as.numeric(pull.pheno(mcross)[, "bw6"]), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

idx <- which(cdata[, "longevity"] >= 365)
cdata <- cdata[idx,]
gtsM <- gtsp[idx,]

minAge <- c(minAge, min(cdata[, "longevity"]))

lods.c <- c()
lm.null <- lm(bw6 ~ sex + site + cohort + 0, data = cdata)
for(marker in colnames(pull.geno(mcross))){
  mp <- gtsM[, grep(marker, colnames(gtsM))]
  lm.alt <- lm(bw6 ~ sex + site + cohort + mp + 0, data = cdata)
  n <- sum(!is.na(lm.alt$resid))
  lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
  lods.c <- c(lods.c, lod)
}
names(lods.c) <- colnames(pull.geno(mcross))


# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

subset <- map[which(map[,1] %in% c(1:19, "X")),]
subset <- cbind(subset, cpos = NA)
gap <- 30000000
chr.start <- c(0)
chr.mids <- c()
cp <- 0
for(chr in c(1:19, "X")){
  cl <- max(as.numeric(subset[which(subset[,1] == chr),2]))
  chr.start <- c(chr.start, cl + cp + gap)
  subset[which(subset[,1] == chr), "cpos"] <- as.numeric(subset[which(subset[,1] == chr), 2]) + cp
  if(chr == "X"){ chr <- 20 }
  cat(chr, " ", cl, " ", cp, " ", gap, " ", chr.start[chr], "\n")
  chr.mids <- c(chr.mids, chr.start[as.numeric(chr)] + cl/2)
  cp = cl + cp + gap
}

plot(c(0, max(chr.start)), y = c(0, 25), t = 'n', ylab = "LOD", xlab = "Chromosome",xaxt="n", las=2, main = "Bodyweight 6mo (≥ 365 days)")
for(x in c(1:19, "X")){
  points(subset[which(subset[,1] == x),"cpos"], lods.c[rownames(subset)[which(subset[,1] == x)]], t = 'l', col = "black",lwd=2)
}
abline(h = 3.65, lty=2, col = "green")
axis(1, at = chr.mids, paste0("", c(1:19, "X")), cex.axis=1.2, las=1)
legend("topleft", c("All", "Males", "Females", "Interaction"), lwd=2, lty=c(1,1,1,3), col = c("black", "blue", "hotpink", "orange"))


genotypes <- pull.geno(fill.geno(mcross, method = "maxmarginal", error.prob=0.01, min.prob = 0.9))
bw6 <- pull.pheno(mcross)[, "bw6"]
sex <- as.factor(pull.pheno(mcross)[, "sex"])
site <- as.factor(pull.pheno(mcross)[, "site"])
cohort <- as.factor(pull.pheno(mcross)[, "cohort"])

lm.null.bw6 <- lm(bw6 ~ sex + site + cohort, data = cdata)
bw6adj <- round(as.numeric(coef(lm.null.bw6)["(Intercept)"]) + residuals(lm.null.bw6), 1)


plot(bw6adj ~ genotypes[,"17_33065600"])


