
library(svglite)
library(pspline)


setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

pos <- "4_55012301"
name <- "vita4a"

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
msequence <- seq(1, max(phe), 15)
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
#curves <- log10(curves)

plot(c(500, 1300), c(0, 1), t = 'n', xlab = "days", ylab = "log10(% surviving)", main = paste0("Survival curve at ", name))
rect(0,0,2000, 3, col = rgb(1,0.8,0.8,0.2))
rect(0,0,2000, -3, col = rgb(0,0.0,0.8,0.2))

curves[is.na(curves)] <- 0


oa.females <- predict(sm.spline(msequence, curves[,2]), msequence, 1)
oa.males <- predict(sm.spline(msequence, curves[,7]), msequence, 1)

op <- par(mfrow=c(2,1))
plot(c(365, 1050), c(0.1, -0.1), t = 'n', xlab = "days", ylab = "Slope - Slope(μ)", main = paste0("Survival curve slope at ", name, " (Female)"))

# Females
#points(msequence, oa, t = 'l', col = "black", lwd=2)
points(msequence, predict(sm.spline(msequence, curves[,3]), msequence, 1) - oa.females, t = 'l', col = "red", lwd=2)
points(msequence, predict(sm.spline(msequence, curves[,4]), msequence, 1) - oa.females, t = 'l', col = "green", lwd=2)
points(msequence, predict(sm.spline(msequence, curves[,5]), msequence, 1) - oa.females, t = 'l', col = "blue", lwd=2)
points(msequence, predict(sm.spline(msequence, curves[,6]), msequence, 1) - oa.females, t = 'l', col = "orange", lwd=2)

plot(c(365, 1050), c(0.1, -0.1), t = 'n', xlab = "days", ylab = "Slope - Slope(μ)", main = paste0("Survival curve slope at ", name, " (Male)"))
# Males
#points(msequence, oa, t = 'l', col = "black", lwd=2)
points(msequence, predict(sm.spline(msequence, curves[,8]), msequence, 1) - oa.males, t = 'l', col = "red", lwd=2)
points(msequence, predict(sm.spline(msequence, curves[,9]), msequence, 1) - oa.males, t = 'l', col = "green", lwd=2)
points(msequence, predict(sm.spline(msequence, curves[,10]), msequence, 1) - oa.males, t = 'l', col = "blue", lwd=2)
points(msequence, predict(sm.spline(msequence, curves[,11]), msequence, 1) - oa.males, t = 'l', col = "orange", lwd=2)


#legend("topright", c("Avg", "C||H", "C||D", "B||H", "B||D"), col = c("black", "red", "green", "blue", "orange"), lwd=2)


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

library(RColorBrewer)
#col.main <- brewer.pal(4, "PiYG")
col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
col.alpha <- add.alpha(col.main, 0.1)

all <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085", "X_150646933")

names(all) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa", "VitaXb")

setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

lods.m.All <- read.table("progressiveMapping_males20D.txt", sep = "\t", check.names=FALSE)
#lods.m.All <- lods.m.All[c(TRUE,FALSE,FALSE), ]
lods.f.All <- read.table("progressiveMapping_females20D.txt", sep = "\t", check.names=FALSE)
#lods.f.All <- lods.f.All[c(TRUE,FALSE,FALSE), ]
lods.c.All <- read.table("progressiveMapping_all20D.txt", sep = "\t", check.names=FALSE)


#### TODO: Update the X to have only 2 haplotypes Males = (C- B-), Females = (CH, BH) (DONE 17/Sept)
#### TODO: Also create a version which starts truncation from 0 days (and not 365 days) (DONE 17/Sept)
#### TODO: Truncation levels for chr 3 (A = T1085 / B = T545) and 4 (T680, Table Max Main effect Table) (DONE 19/Sept)
#### TODO: Add Sex x G plot (Delta males versus females) (Done 25/Sept)
#### TODO: Add text about inversions / recombinations in M&M

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__bioRxiv_All_Key_Files/11_FiguresDanny")

for(ii in 1:length(all)){
  pos <- all[ii]
  name = names(all)[ii]

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

  rownames(remaining) <- paste0("day", msequence)
  rownames(errors) <- paste0("day", msequence)

  pdf(paste0(name, ".eff.new.pdf"), width = 36, height = 12)
  op <- par(mfrow = c(1,3))
  op <- par(cex = 2)
  par(mar = c(5, 5, 4, 3))
  plot(c(20, 1100), c(-40, 35), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
         ylab = "Actuarial Effect Size [d]", main = paste0(name," Combined"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
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

  points(msequence, remaining[,4], t = 'l', col = col.main[3], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,4] + errors[,3], rev(remaining[,4] - errors[,3])), col = col.alpha[3], border = NA)

  points(msequence, remaining[,5], t = 'l', col = col.main[4], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,5] + errors[,4], rev(remaining[,5] - errors[,4])), col = col.alpha[4], border = NA)

  #rect(0, -50, 1200, -30, col = "white", border=NA)
  abline(h = -40 + 4.65, lty = 1, col = "green")
  abline(h = -40 + 3.95, lty = 1, col = "orange")
  abline(h = -40 + 3, lty = 1, col = "red")
  points(seq(20, 1100, 15), lods.c.All[, all[ii]] - 40, t = "l", lwd=2)
  axis(4, at = seq(-40, -32, 2), seq(0, 8, 2), las = 2,  cex.axis =1.1)

  #abline(v=c(935, 1055))
  if(grepl("X", name)){
    legend("top", c("CH", "C-", "BH", "B-"), col = col.main, lwd=4, bg = "white", ncol=4, bty = "n")
  }else{
    legend("top", c("CH", "CD", "BH", "BD"), col = col.main, lwd=4, bg = "white", ncol=4, bty = "n")
  }
  box()

  mtext('LOD', side=4, line=0.7, at=-26, cex = 2.8)


  # Females
  plot(c(20, 1100), c(-40, 35), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
        ylab = "", main = paste0(name," Females"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2, cex.axis = 1.4)
  #abline(v = seq(400, 1100, 50), lty=2, col="lightgray")
  #abline(h = seq(-50, 50, 10), lty=2, col="lightgray")

  points(msequence, remaining[,7], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,7] + errors[,5], rev(remaining[,7] - errors[,5])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,8], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,8] + errors[,6], rev(remaining[,8] - errors[,6])), col = col.alpha[2], border = NA)

  points(msequence, remaining[,9], t = 'l', col = col.main[3], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,9] + errors[,7], rev(remaining[,9] - errors[,7])), col = col.alpha[3], border = NA)

  points(msequence, remaining[,10], t = 'l', col = col.main[4], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,10] + errors[,8], rev(remaining[,10] - errors[,8])), col = col.alpha[4], border = NA)
  #abline(v=1040)

  #rect(0, -50, 1200, -30, col = "white", border=NA)
  abline(h = -40 + 4.65, lty = 1, col = "green")
  abline(h = -40 + 3.95, lty = 1, col = "orange")
  abline(h = -40 + 3, lty = 1, col = "red")
  points(seq(20, 1100, 15), lods.f.All[, all[ii]] - 40, t = "l", lwd=2)
  axis(4, at = seq(-40, -32, 2), seq(0, 8, 2), las = 2,  cex.axis =1.1)
  if(grepl("X", name)){
    legend("top", c("CH", "CD", "BH", "BD")[c(1,3)], col = col.main[c(1,3)], lwd=4, bg = "white", ncol=4, bty = "n")
  }else{
    legend("top", c("CH", "CD", "BH", "BD"), col = col.main, lwd=4, bg = "white", ncol=4, bty = "n")
  }
  box()

  mtext('LOD', side=4, line=0.7, at=-26, cex = 2.8)

  # Males
  plot(c(20, 1100), c(-40, 35), t = 'n', xlab = "Truncation Age [d]", yaxs = "i",
         ylab = "", main = paste0(name," Males"), xaxs = "i", cex.main = 1.4, cex.lab = 1.4, las=2, xaxt="n", yaxt="n")
  axis(1, at = c(20, 200, 400, 600, 800, 1000), c(20, 200, 400, 600, 800, 1000),cex.axis = 1.4)
  axis(1, at = c(100, 300, 500, 700, 900, 1100), c("", "","", "", "", ""))
  aa <- seq(-50, 50, 5)
  aa[which(1:length(aa) %% 2 == 0)] <- ""
  axis(2, at = seq(-50, 50, 5), aa, las=2,cex.axis = 1.4)
  #abline(v = seq(400, 1100, 50), lty=2, col="lightgray")
  #abline(h = seq(-50, 50, 10), lty=2, col="lightgray")

  points(msequence, remaining[,12], t = 'l', col = col.main[1], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,12] + errors[,9], rev(remaining[,12] - errors[,9])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,13], t = 'l', col = col.main[2], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,13] + errors[,10], rev(remaining[,13] - errors[,10])), col = col.alpha[2], border = NA)

  points(msequence, remaining[,14], t = 'l', col = col.main[3], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,14] + errors[,11], rev(remaining[,14] - errors[,11])), col = col.alpha[3], border = NA)

  points(msequence, remaining[,15], t = 'l', col = col.main[4], lwd=3)
  polygon(c(msequence, rev(msequence)), c(remaining[,15] + errors[,12], rev(remaining[,15] - errors[,12])), col = col.alpha[4], border = NA)

  #rect(0, -50, 1200, -30, col = "white", border=NA)
  abline(h = -40 + 4.65, lty = 1, col = "green")
  abline(h = -40 + 3.95, lty = 1, col = "orange")
  abline(h = -40 + 3, lty = 1, col = "red")
  points(seq(20, 1100, 15), lods.m.All[, all[ii]] - 40, t = "l", lwd=2)
  axis(4, at = seq(-40, -32, 2), seq(0, 8, 2), las = 2,  cex.axis =1.1)
  if(grepl("X", name)){
    legend("top", c("CH", "C-", "BH", "B-")[c(2,4)], col = col.main[c(2,4)], lwd=4, bg = "white", ncol=4, bty = "n")
  }else{
    legend("top", c("CH", "CD", "BH", "BD"), col = col.main, lwd=4, bg = "white", ncol=4, bty = "n")
  }
  box()

  mtext('LOD', side=4, line=0.7, at=-26, cex = 2.8)

  dev.off()
}


### PLot 3 x 3

all <- all[c("Vita1b", "Vita9c", "Vita15a")]
svglite(paste0("paper.eff.svg"), width = (1200 * (600/72)) / 280, height = (1200 * (600/72)) / 280)
nf <- layout(matrix(c(1,2,3,4,
                      5,6,7,8,
                      9,10,11,12,
                      13,14,15,16), ncol=4, byrow=TRUE), heights=c(1, 12, 12, 12), widths=c(2,5,5,5))

par(cex=1.4)

  par(mar = c(2, 2, 2, 2))
  plot.new()
  text(0.5,0.5,"", cex=1.4, font=2, srt = 90)
  legend("center", c("C||H", "C||D", "B||H", "B||D"), col = col.main, lwd=2, bg = "white", cex=0.7, ncol=2, bty = 'n')

  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5,0.5,"Combined", cex=1.4, font=2)

  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5,0.5,"Females", cex=1.4, font=2)

  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5,0.5,"Males", cex=1.4, font=2)


for(ii in 1:1){ #length(all)){
  pos <- all[ii]
  name = names(all)[ii]

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

  par(mar = c(2, 2, 2, 2))
  plot.new()
  text(0.5,0.5,substitute(italic(x), list(x = name)), cex=1.4, font=2, srt = 90)


  par(mar = c(2.1, 2.5, 1.1, 0.5))
  plot(c(365, 1050), c(-40, 25), t = 'n', xlab = "", ylab = "Allele effect (days)", main = "", xaxs = "i", las=1, xaxt="n")
  axis(1, at = c(400, 600, 800, 1000), c(400, 600, 800, 1000))
  abline(v = seq(400, 1100, 50), lty=2, col="lightgray")
  abline(h = seq(-50, 50, 10), lty=2, col="lightgray")
  text(x=240, y=0, 'Allele effect (days)', xpd=NA, srt = 90)

  # Combined
  points(msequence, remaining[,2], t = 'l', col = col.main[1], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,2] + errors[,1], rev(remaining[,2] - errors[,1])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,3], t = 'l', col = col.main[2], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,3] + errors[,2], rev(remaining[,3] - errors[,2])), col = col.alpha[2], border = NA)

  points(msequence, remaining[,4], t = 'l', col = col.main[3], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,4] + errors[,3], rev(remaining[,4] - errors[,3])), col = col.alpha[3], border = NA)

  points(msequence, remaining[,5], t = 'l', col = col.main[4], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,5] + errors[,4], rev(remaining[,5] - errors[,4])), col = col.alpha[4], border = NA)

  #abline(v=c(935, 1055))
  #legend("topleft", c("C||H", "C||D", "B||H", "B||D"), col = col.main, lwd=2, bg = "white", cex=0.7, ncol=2)
  box()

  par(mar = c(2.1, 2.5, 1.1, 0.5))

  # Females
  plot(c(365, 1050), c(-40, 25), t = 'n', xlab = "", ylab = "Allele effect (days)", main = "", xaxs = "i", las=1, xaxt="n")
  axis(1, at = c(400, 600, 800, 1000), c(400, 600, 800, 1000))
  abline(v = seq(400, 1100, 50), lty=2, col="lightgray")
  abline(h = seq(-50, 50, 10), lty=2, col="lightgray")

  points(msequence, remaining[,7], t = 'l', col = col.main[1], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,7] + errors[,5], rev(remaining[,7] - errors[,5])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,8], t = 'l', col = col.main[2], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,8] + errors[,6], rev(remaining[,8] - errors[,6])), col = col.alpha[2], border = NA)

  points(msequence, remaining[,9], t = 'l', col = col.main[3], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,9] + errors[,7], rev(remaining[,9] - errors[,7])), col = col.alpha[3], border = NA)

  points(msequence, remaining[,10], t = 'l', col = col.main[4], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,10] + errors[,8], rev(remaining[,10] - errors[,8])), col = col.alpha[4], border = NA)
  #abline(v=1040)
  #legend("topleft", c("C||H", "C||D", "B||H", "B||D"), col = col.main, lwd=2, bg = "white", cex=0.7, ncol=2)
  box()

  par(mar = c(2.1, 2.5, 1.1, 0.5))
  # Males
  plot(c(365, 1050), c(-40, 25), t = 'n', xlab = "", ylab = "Allele effect (days)", main = "", xaxs = "i", las=1, xaxt="n")
  axis(1, at = c(400, 600, 800, 1000), c(400, 600, 800, 1000))
  abline(v = seq(400, 1100, 50), lty=2, col="lightgray")
  abline(h = seq(-50, 50, 10), lty=2, col="lightgray")

  points(msequence, remaining[,12], t = 'l', col = col.main[1], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,12] + errors[,9], rev(remaining[,12] - errors[,9])), col = col.alpha[1], border = NA)

  points(msequence, remaining[,13], t = 'l', col = col.main[2], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,13] + errors[,10], rev(remaining[,13] - errors[,10])), col = col.alpha[2], border = NA)

  points(msequence, remaining[,14], t = 'l', col = col.main[3], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,14] + errors[,11], rev(remaining[,14] - errors[,11])), col = col.alpha[3], border = NA)

  points(msequence, remaining[,15], t = 'l', col = col.main[4], lwd=2)
  polygon(c(msequence, rev(msequence)), c(remaining[,15] + errors[,12], rev(remaining[,15] - errors[,12])), col = col.alpha[4], border = NA)
  box()

  #legend("topleft", c("C||H", "C||D", "B||H", "B||D"), col = col.main, lwd=2, bg = "white", cex=0.7, ncol=2)
}
dev.off()



