setwd("D:/Ddrive/Github/UM-HET3")
source("adjustXprobs.R")
setwd("C:/Github/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)
nxo <- countXO(mcross)

cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                    bw6 = pull.pheno(mcross)[, "bw6"], 
                    bw12 = pull.pheno(mcross)[, "bw12"], 
                    bw18 = pull.pheno(mcross)[, "bw18"], 
                    bw24 = pull.pheno(mcross)[, "bw24"], 
                    nxo = as.numeric(nxo), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
idx <- which(cdata[, "longevity"] > 365)
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

# Beta, std.err and pvalues of covariate effects
summary(lm(longevity ~ sex + 0, data = cdata))
anova(lm(longevity ~ sex, data = cdata))

summary(lm(longevity ~ site + 0, data = cdata))
anova(lm(longevity ~ site, data = cdata))

summary(lm(longevity ~ cohort + 0, data = cdata))
anova(lm(longevity ~ cohort, data = cdata))

summary(lm(longevity ~ treatment + 0, data = cdata))
anova(lm(longevity ~ treatment, data = cdata))

cor(cdata[,1:5], use="pair", method = "spearman")

## Combined longevity > 356 days
log10.I <- matrix(NA, ncol(pull.geno(mcross)), ncol(pull.geno(mcross)), dimnames=list(colnames(pull.geno(mcross)), colnames(pull.geno(mcross))))
for(m1 in colnames(pull.geno(mcross))){
  mp1 <- gtsp[, grep(m1, colnames(gtsp))]
  for(m2 in colnames(pull.geno(mcross))){
    mp2 <- gtsp[, grep(m2, colnames(gtsp))]
    lm.null <- lm(longevity ~ sex + site + cohort + treatment + mp1 + mp2 + 0, data = cdata)
    lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp1 * mp2 + 0, data = cdata)
    log10.I[m1,m2] <- -log10(anova(lm.alt, lm.null)[2,"Pr(>F)"])
  }
}
write.table(log10.I, "GxG_LODs.txt", sep = "\t")

## Combined longevity > 356 days
lods.c <- c()
lm.null <- lm(longevity ~ sex + site + cohort + treatment + 0, data = cdata)
for(marker in colnames(pull.geno(mcross))){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp + 0, data = cdata)
  n <- sum(!is.na(lm.alt$resid))
  lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
  lods.c <- c(lods.c, lod)
}
names(lods.c) <- colnames(pull.geno(mcross))


## Females
females <- which(cdata[, "sex"] == 0)

lods.f <- c()
lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = cdata[females, ])
for(marker in colnames(pull.geno(mcross))){
  mp <- gtsp[females, grep(marker, colnames(gtsp))]
  lm.alt <- lm(longevity ~ site + cohort + treatment + mp + 0, data = cdata[females, ])
  n <- sum(!is.na(lm.alt$resid))
  lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
  lods.f <- c(lods.f, lod)
}
names(lods.f) <- colnames(pull.geno(mcross))

## Males
males <- which(cdata[, "sex"] == 1)

lods.m <- c()
lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = cdata[males, ])
for(marker in colnames(pull.geno(mcross))){
  mp <- gtsp[males, grep(marker, colnames(gtsp))]
  lm.alt <- lm(longevity ~ site + cohort + treatment + mp + 0, data = cdata[males, ])
  n <- sum(!is.na(lm.alt$resid))
  lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
  lods.m <- c(lods.m, lod)
}
names(lods.m) <- colnames(pull.geno(mcross))

lods.all <- cbind(lods.c, lods.f, lods.m)
write.table(lods.all, file = "lods.all.txt")
significant <- unlist(lapply(strsplit(names(which(apply(lods.all,1, function(x){any(x > 4.25)}))), "_"), "[",1))

# Plot the QTL profile
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))
chr.l <- c()
for(chr in unique(chrs)){ chr.l <- c(chr.l, max(positions[which(chrs==chr)])); }
names(chr.l) <- unique(chrs)

plot(c(-1000000, 5000000 + sum(chr.l) + 19 * 25000000), c(0,8), main = "QTL profiles for longevity", t = 'n', xaxt='n', xlab="Chromosome", ylab="LOD", xaxs="i", yaxs="i", las=2)
#abline(h = seq(0,10,2), lty=2)
s <- 0
h <- c()
for(chr in unique(chrs)){
  rect(s, 0, s + chr.l[chr] + 5000000, 12, col = rgb(0.9,0.9,0.9,0.5),border=NA)
  points(positions[chrs==chr] + s, lods.c[chrs == chr], t = 'l', lwd=3, col="black")
  points(positions[chrs==chr] + s, lods.f[chrs == chr], t = 'l', lwd=2, col="pink")
  points(positions[chrs==chr] + s, lods.m[chrs == chr], t = 'l', lwd=2, col="blue")
  h <- c(h, s + (chr.l[chr]/2) + 25000000)
  s <- s + chr.l[chr] + 25000000
}
axis(1, at = h, unique(chrs))
abline(h = c(4.25, 4.95), col = c("orange", "green"), lty=2)
legend("topright", c("Males and females combined", "Only females", "Only males"), lwd=c(3,2,2), col=c("black", "pink", "blue"), bg="white")
legend("topleft", c("p = 0.05", "p = 0.01"), lwd=1, lty=2, col=c("orange", "green"), bg="white")
box()

#QTL main effects on adjusted longevity
getEffect <- function(sdata, gtsprob, marker = "1_24042124", model = "longevity ~ sex + site + cohort + treatment"){
  rownames(sdata) <- 1:nrow(sdata)
  rownames(gtsprob) <- 1:nrow(gtsprob)
  mp <- gtsprob[, grep(marker, colnames(gtsprob))]
  gts <- unlist(lapply(lapply(lapply(apply(mp,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  mlm <- lm(as.formula(model), data = sdata)
  pheAdj <- rep(NA, nrow(sdata))
  adj <- residuals(mlm) + mean(sdata[, "longevity"])
  pheAdj[as.numeric(names(adj))] <- adj
  OAmean <- mean(pheAdj[which(!is.na(gts))])
  means <- c(mean(pheAdj[which(gts == "AC")]),mean(pheAdj[which(gts == "AD")]),mean(pheAdj[which(gts == "BC")]),mean(pheAdj[which(gts == "BD")])) - OAmean
  std <- function(x) sd(x)/sqrt(length(x))
  stderrs <- c(std(pheAdj[which(gts == "AC")]),std(pheAdj[which(gts == "AD")]),std(pheAdj[which(gts == "BC")]),std(pheAdj[which(gts == "BD")]))
  paste0(length(which(!is.na(gts))), "\t", round(OAmean,0), "\t", paste0(round(means,0), " ± ", round(stderrs,2), collapse="\t"))
}

write.table(
rbind(getEffect(cdata, gtsp, "1_3010274"),
      getEffect(cdata, gtsp, "4_88722326"),
      getEffect(cdata, gtsp, "X_36008085"))
      , file = "All_M+F_Eff.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(
rbind(getEffect(cdata[females,], gtsp[females,], "2_144490178", "longevity ~ site + cohort + treatment"),
      getEffect(cdata[females,], gtsp[females,], "9_34932404", "longevity ~ site + cohort + treatment"))
      , file = "All_F_Eff.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(
rbind(getEffect(cdata[males,], gtsp[males,], "6_134870385", "longevity ~ site + cohort + treatment"))
      , file = "All_M_Eff.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


mdata <- cdata[males,]
gtspm <- gtsp[males,]
idxs <- sort(mdata[, "longevity"], index.return=TRUE)$ix
getEffect(mdata[idxs[(1+ (.42 * length(idxs))): length(idxs)],], gtspm[idxs[(1+ (.42 * length(idxs))): length(idxs)],], "1_3010274", "longevity ~ site + cohort + treatment") #58%
min(mdata[idxs[(1+ (.42 * length(idxs))): length(idxs)],"longevity"])

getEffect(mdata[idxs[(1+ (.20 * length(idxs))): length(idxs)],], gtspm[idxs[(1+ (.20 * length(idxs))): length(idxs)],], "4_52524395", "longevity ~ site + cohort + treatment") #80%
getEffect(mdata[idxs[(1+ (.38 * length(idxs))): length(idxs)],], gtspm[idxs[(1+ (.38 * length(idxs))): length(idxs)],], "9_124056586", "longevity ~ site + cohort + treatment") #62%
getEffect(mdata[idxs[(1+ (.74 * length(idxs))): length(idxs)],], gtspm[idxs[(1+ (.74 * length(idxs))): length(idxs)],], "10_72780332", "longevity ~ site + cohort + treatment") #26%
getEffect(mdata[idxs[(1+ (.42 * length(idxs))): length(idxs)],], gtspm[idxs[(1+ (.42 * length(idxs))): length(idxs)],], "15_55481391", "longevity ~ site + cohort + treatment") #58%
getEffect(mdata[idxs[(1+ (.26 * length(idxs))): length(idxs)],], gtspm[idxs[(1+ (.26 * length(idxs))): length(idxs)],], "17_34460077", "longevity ~ site + cohort + treatment") #74%

fdata <- cdata[females,]
gtspf <- gtsp[females,]
idxs <- sort(fdata[, "longevity"], index.return=TRUE)$ix
getEffect(fdata[idxs[(1+ (.22 * length(idxs))): length(idxs)],], gtspf[idxs[(1+ (.22 * length(idxs))): length(idxs)],], "1_24042124", "longevity ~ site + cohort + treatment") #78%
min(fdata[idxs[(1+ (.42 * length(idxs))): length(idxs)],"longevity"])

getEffect(fdata[idxs[(1+ (.02 * length(idxs))): length(idxs)],], gtspf[idxs[(1+ (.02 * length(idxs))): length(idxs)],], "2_148442635", "longevity ~ site + cohort + treatment") #98%
getEffect(fdata[idxs[(1+ (.08 * length(idxs))): length(idxs)],], gtspf[idxs[(1+ (.08 * length(idxs))): length(idxs)],], "9_29939029", "longevity ~ site + cohort + treatment") #92%
getEffect(fdata[idxs[(1+ (.16 * length(idxs))): length(idxs)],], gtspf[idxs[(1+ (.16 * length(idxs))): length(idxs)],], "X_156343080", "longevity ~ site + cohort + treatment") #84%

idxs <- sort(cdata[, "longevity"], index.return=TRUE)$ix
getEffect(cdata[idxs[(1+ (.20 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.20 * length(idxs))): length(idxs)],], "1_24042124", "longevity ~ sex + site + cohort + treatment") #80%
getEffect(cdata[idxs[(1+ (.34 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.34 * length(idxs))): length(idxs)],], "2_112712327", "longevity ~ sex + site + cohort + treatment") #66%
getEffect(cdata[idxs[(1+ (.14 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.14 * length(idxs))): length(idxs)],], "4_55012301", "longevity ~ sex + site + cohort + treatment") #86%
getEffect(cdata[idxs[(1+ (.04 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.04 * length(idxs))): length(idxs)],], "6_107382038", "longevity ~ sex + site + cohort + treatment") #96%
getEffect(cdata[idxs[(1+ (.72 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.72 * length(idxs))): length(idxs)],], "9_122888918", "longevity ~ sex + site + cohort + treatment") #28%
getEffect(cdata[idxs[(1+ (.70 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.70 * length(idxs))): length(idxs)],], "10_72780332", "longevity ~ sex + site + cohort + treatment") #30%
getEffect(cdata[idxs[(1+ (.14 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.14 * length(idxs))): length(idxs)],], "12_112855820", "longevity ~ sex + site + cohort + treatment") #86%
getEffect(cdata[idxs[(1+ (.46 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.46 * length(idxs))): length(idxs)],], "14_101437457", "longevity ~ sex + site + cohort + treatment") #54%
getEffect(cdata[idxs[(1+ (.54 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.54 * length(idxs))): length(idxs)],], "15_74248242", "longevity ~ sex + site + cohort + treatment") #46%
getEffect(cdata[idxs[(1+ (.1 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.1 * length(idxs))): length(idxs)],], "17_32883804", "longevity ~ sex + site + cohort + treatment") #90%
getEffect(cdata[idxs[(1+ (.24 * length(idxs))): length(idxs)],], gtsp[idxs[(1+ (.24 * length(idxs))): length(idxs)],], "X_156343080", "longevity ~ sex + site + cohort + treatment") #76%





