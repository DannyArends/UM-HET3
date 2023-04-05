setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

pos <- "9_104091597"
name = "vita9c"

mm <- pull.geno(fill.geno(mcross))
phe <- pull.pheno(mcross)[, "longevity"]

sex <- pull.pheno(mcross)[, "sex"]
site <- as.factor(pull.pheno(mcross)[, "site"])
cohort <- as.factor(pull.pheno(mcross)[, "cohort"])
treatment <- as.factor(pull.pheno(mcross)[, "treatment"])

marker <- mm[,pos]

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

plot(c(500, 1300), c(-2, 2), t = 'n', xlab = "days", ylab = "log10(% surviving)", main = paste0("Survival curve at ", name))
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
  means <- c(mean(pheAdj[which(gts == "AC")]),mean(pheAdj[which(gts == "AD")]),mean(pheAdj[which(gts == "BC")]),mean(pheAdj[which(gts == "BD")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "AC")]),std(pheAdj[which(gts == "AD")]),std(pheAdj[which(gts == "BC")]),std(pheAdj[which(gts == "BD")]))
  #paste0(length(which(!is.na(gts))), "/",nrow(cdata), "\t", round(OAmean,0), "\t", paste0(round(means,0), " ± ", round(stderrs,2), collapse="\t"))
  return(list(
          c( round(OAmean,0), round(means,0)),
          c(round(stderrs,2))
        ))
}

remaining <- c()
errors <- c()
msequence <- seq(365, 1100, 15)
for(d in msequence){
  combined <- getEffect(mcross, gtsp, marker = pos, timepoint = d)
  female <- getEffect(mcross, gtsp, marker = pos, timepoint = d, sex = 0, model = "longevity ~ site + cohort + treatment")
  male <- getEffect(mcross, gtsp, marker = pos, timepoint = d, sex = 1, model = "longevity ~ site + cohort + treatment")
  remaining <- rbind(remaining, c(combined[[1]], female[[1]], male[[1]]))
  errors <- rbind(errors, c(combined[[2]], female[[2]], male[[2]]))
}

rownames(remaining) <- paste0("day", msequence)
rownames(errors) <- paste0("day", msequence)

op <- par(mfrow = c(1,3))
plot(c(365, 1100), c(-40, 40), t = 'n', xlab = "days", ylab = "Allele effect (days)", main = paste0("Allele effect at ",name," (combined)"))

# Combined
points(msequence, remaining[,2], t = 'l', col = "red", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,2] + errors[,1], rev(remaining[,2] - errors[,1])), col = rgb(1,0,0,0.2), border = NA)

points(msequence, remaining[,3], t = 'l', col = "green", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,3] + errors[,2], rev(remaining[,3] - errors[,2])), col = rgb(0,1,0,0.2), border = NA)

points(msequence, remaining[,4], t = 'l', col = "blue", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,4] + errors[,3], rev(remaining[,4] - errors[,3])), col = rgb(0,0,1,0.2), border = NA)

points(msequence, remaining[,5], t = 'l', col = "orange", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,5] + errors[,4], rev(remaining[,5] - errors[,4])), col = rgb(1,.647, 0, 0.2), border = NA)

#abline(v=c(935, 1055))
legend("topleft", c("C||H", "C||D", "B||H", "B||D"), col = c("red", "green", "blue", "orange"), lwd=2, bg = "white")

plot(c(365, 1100), c(-40, 40), t = 'n', xlab = "days", ylab = "Allele effect (days)", main = paste0("Allele effect at ",name," (females)"))
# Females
points(msequence, remaining[,7], t = 'l', col = "red", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,7] + errors[,5], rev(remaining[,7] - errors[,5])), col = rgb(1,0,0,0.2), border = NA)

points(msequence, remaining[,8], t = 'l', col = "green", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,8] + errors[,6], rev(remaining[,8] - errors[,6])), col = rgb(0,1,0,0.2), border = NA)

points(msequence, remaining[,9], t = 'l', col = "blue", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,9] + errors[,7], rev(remaining[,9] - errors[,7])), col = rgb(0,0,1,0.2), border = NA)

points(msequence, remaining[,10], t = 'l', col = "orange", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,10] + errors[,8], rev(remaining[,10] - errors[,8])), col = rgb(1,.647, 0, 0.2), border = NA)
#abline(v=1040)
legend("topleft", c("C||H", "C||D", "B||H", "B||D"), col = c("red", "green", "blue", "orange"), lwd=2, bg = "white")

# Males
plot(c(365, 1100), c(-40, 40), t = 'n', xlab = "days", ylab = "Allele effect (days)", main = paste0("Allele effect at ",name," (males)"))
points(msequence, remaining[,12], t = 'l', col = "red", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,12] + errors[,9], rev(remaining[,12] - errors[,9])), col = rgb(1,0,0,0.2), border = NA)

points(msequence, remaining[,13], t = 'l', col = "green", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,13] + errors[,10], rev(remaining[,13] - errors[,10])), col = rgb(0,1,0,0.2), border = NA)

points(msequence, remaining[,14], t = 'l', col = "blue", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,14] + errors[,11], rev(remaining[,14] - errors[,11])), col = rgb(0,0,1,0.2), border = NA)

points(msequence, remaining[,15], t = 'l', col = "orange", lwd=2)
polygon(c(msequence, rev(msequence)), c(remaining[,15] + errors[,12], rev(remaining[,15] - errors[,12])), col = rgb(1,.647, 0, 0.2), border = NA)

legend("topleft", c("C||H", "C||D", "B||H", "B||D"), col = c("red", "green", "blue", "orange"), lwd=2, bg = "white")


