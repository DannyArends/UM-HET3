#
# ctl_sex.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Visualization (in CTL style) of Male/Female correlation differences using 6 month bodyweight (think: sex = genetic marker)
#

library(qtl)

source("ActuarialMapping/adjustXprobs.R")
mcross <- read.cross(format="csvr", file="DataSet/um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
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
idx <- which(cdata[, "longevity"] >= (365/2) & !is.na(cdata[, "bw6"]))
cdata <- cdata[idx,]
gtsp <- gtsp[idx,]

lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjLongevity"] <- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)

lm.null.bw6 <- lm(bw6 ~ sex + site + cohort + treatment, data = cdata)
cdata[, "adjBw6"] <- round(as.numeric(coef(lm.null.bw6)["(Intercept)"]) + residuals(lm.null.bw6), 2)

lm(cdata[which(cdata[, "sex"] == 0), "adjLongevity"] ~ cdata[which(cdata[, "sex"] == 0), "adjBw6"]+1);

corM <- c()
allN <- c()
confL <- c()
confU <- c()
for(x in seq(185, 1100, 15)){
  male <- which(cdata[, "sex"] %in% 0 & cdata[, "adjLongevity"] > x)
  fema <- which(cdata[, "sex"] %in% 1 & cdata[, "adjLongevity"] > x)
  cMale <- NA; cFema <- NA;

  if(length(male) > 100){
    cMale <- cor(cdata[male, "adjLongevity"], cdata[male, "adjBw6"], use = "pair", method = "spearman");
    confMale <- cor.test(cdata[male, "adjLongevity"], cdata[male, "adjBw6"], use = "pair", method = "pearson", conf.level = 0.50);
  }
  if(length(fema) > 100){
    cFema <- cor(cdata[fema, "adjLongevity"], cdata[fema, "adjBw6"], use = "pair", method = "spearman");
    confFema <- cor.test(cdata[fema, "adjLongevity"], cdata[fema, "adjBw6"], use = "pair", method = "pearson", conf.level = 0.50);
  }

  corM <- rbind(corM, c(cMale, cFema))

  confL <- rbind(confL, c(round(as.numeric(unlist(confMale)["conf.int1"]),2),
                          round(as.numeric(unlist(confFema)["conf.int1"]),2)))
  confU <- rbind(confU, c(round(as.numeric(unlist(confMale)["conf.int2"]),2),
                          round(as.numeric(unlist(confFema)["conf.int2"]),2)))
  # Adjust Conf
  for(x in 1:nrow(corM)){
    for(y in 1:ncol(corM)){
      mid <- corM[x, y]
      adj <- (confU[x,y] - confL[x, y]) / 2

      confU[x,y] <- mid + adj
      if(is.na(confU[x,y])) confU[x,y] <- 0
      confL[x,y] <- mid - adj
      if(is.na(confL[x,y])) confL[x,y] <- 0
    }
  }
}
col.main <- c("#FF3333", "#00AEEF")
add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
col.alpha2 <- add.alpha(col.main, 0.1)

plot(c(185, 1100), c(-0.5, 0.1), t = "n", xlab = "Truncation age (days)", ylab = "Correlation BW185 to Tage", 
     main = "Correlation Male/Female 6months BW",yaxt = "n", yaxs = "i")

points(seq(185, 1100, 15), corM[,1], t = "l", col = col.main[1], lwd = 2)
polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,1], rev(confU[,1])),col = col.alpha2[1], border=NA)

points(seq(185, 1100, 15), corM[,2], t = "l", col = col.main[2], lwd = 2)
polygon(c(seq(185, 1100, 15), rev(seq(185, 1100, 15))), c(confL[,2], rev(confU[,2])),col = col.alpha2[2], border=NA)

axis(2, at = seq(-.5, .1, .1), las=2)
axis(2, at = seq(-.5, .1, .05), rep("", length(seq(-.5, .1, .05))), las=2)
axis(1, at = seq(100, 1100, 100), rep("",length(seq(100, 1100, 100))))

legend("topleft", c("Female", "Male"), col = col.main, lwd=4, bg = "white", ncol=4, bty = "n")

