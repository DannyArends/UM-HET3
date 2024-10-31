library(RColorBrewer)
library(svglite)
library(vioplot)

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3")

### Interaction plot
setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

all <- c("1_3010272", "1_24042124", "1_120474787", "2_89844287", "2_112712327", "2_139956785","3_83838529", "3_92135706", "4_52524395",
         "5_67573068", "6_107382038", "6_140010438", "9_29939029", "9_104091597", "9_124056586", "10_72780332", "11_5628810", "11_82178599",
         "12_112855820", "13_89689878", "14_101437457", "15_74248242", "17_32883804", "18_60822951", "X_36008085", "X_156343080")

names(all) <- c("Vita1a","Vita1b","Vita1c","Vita2a","Vita2b","Vita2c","Vita3a","Vita3b","Vita4a","Vita5a","Vita6a","Vita6b",
                "Vita9a","Vita9b","Vita9c","Vita10a","Vita11a","Vita11b","Vita12a","Vita13a","Vita14a","Vita15a",
                "Vita17a","Vita18a","VitaXa","VitaXb")

for(m in 1:length(all)){

  marker1 <- all[m]

  mp1 <- gtsp[, grep(marker1, colnames(gtsp))]
  gts1 <- unlist(lapply(lapply(lapply(apply(mp1,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]),
                      gts1 = gts1)

  idx <- which(cdata[, "longevity"] >= 365)
  cdata <- cdata[idx,]

  mlm <- lm("longevity ~ sex + site + cohort + treatment", data = cdata)

  pheAdj <- rep(NA, nrow(cdata))
  names(pheAdj) <-  rownames(cdata)

  adj <- residuals(mlm) + mean(cdata[, "longevity"])
  pheAdj[names(adj)] <- round(adj)

  cdata <- cbind(pheAdj, cdata)
  cdata <- cdata[which(!is.na(cdata[, "gts1"]) & !is.na(cdata[, "sex"])),]
  std <- function(x) sd(x)/sqrt(length(x))
  col.main <- c("#DB7093", "#004BAD", "#B16BE6", "#F02D27")
  add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
  col.alpha <- add.alpha(col.main, 0.1)

  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_BioRxiv_All_Key_Files/03_Figure_3_GxS")

  cdata <- cdata[which(cdata[, "pheAdj"] > 365),]
  pdf(paste0(names(all)[m], ".GxS.pdf"), width = 12, height = 12)
  op <- par(cex = 2)
  plot(c(1,4), c(mean(cdata[, "pheAdj"])-40, mean(cdata[, "pheAdj"])+40), t = "n", 
       main = paste0(names(all)[m], " x Sex interaction (>=",365," days)"), 
       xlab = paste0("Haplotype ",names(all)[m]), ylab = "Adjusted Longevity (days)", xaxt="n")
  off <- -0.15
  j <- 1
  for(v1 in c(0, 1)){
    i <- 1
    vals <- c()
    errs <- c()
    for(v2 in c("AC", "AD", "BC", "BD")){
      values <- cdata[which(cdata[, "gts1"] == v2 & cdata[, "sex"] == v1), "pheAdj"]
      #vioplot(values, at = i + off, add = TRUE, wex=0.1, col = col.main[j])
      vals <- c(vals, mean(values))
      errs <- c(errs, std(values))
      i <- i + 1
    }
    points(vals, col=col.main[j], t = "b",pch=20)
    polygon(c(1:4, 4:1), c(vals + errs, rev(vals - errs)), col = col.alpha[j], border = NA)
    off <- off + 0.1
    j <- j + 1
  }
  axis(1, at = 1:4, c("CH", "CD", "BH", "BD"))
  legend("topleft", c("Females", "Males"), col = col.main, lwd=1, pch = 20, title = "Sex")
  dev.off()
}

