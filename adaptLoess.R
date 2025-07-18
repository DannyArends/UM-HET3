# TODO: Follow up with Arthus for the supplemental files - Ask again
# TODO: Loess smoother
# TODO: Output a folder of these for all Vita / haplotype
# TODO: M&M Interaction scan 
# TODO: Mortality figure for KM (Figure 1b) Males versus Females whole population


all <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085")

names(all) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa")


setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)

gtsp <- pull.genoprob(mcross)

col.main <- c("#00cdcd", "#9400d3", "#ff1493", "#deb887")
ss <- c("Male", "Female")
names(ss) <- c("1", "0")
for(mname in names(all)){
  marker <- all[mname]

  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  areA <- grep("A", gts)
  areB <- grep("B", gts)
  areC <- grep("C", gts)
  areD <- grep("D", gts)

  gtsM <- rep(NA, length(gts))
  gtsM[areA] <- "C" ## BALB
  gtsM[areB] <- "B" ## B6

  gtsP <- rep(NA, length(gts))
  gtsP[areC] <- "H" ## C3H
  gtsP[areD] <- "D" ## DBA

  phe <- cbind(lifespan = pull.pheno(mcross)[, "longevity"], sex = pull.pheno(mcross)[, "sex"], gtsM, gtsP)

  for(sex in c("0", "1")) {
    phe.m <- phe[which(phe[, "sex"] == sex & !is.na(phe[, "gtsM"])),]
    phe.p <- phe[which(phe[, "sex"] == sex & !is.na(phe[, "gtsP"])),]

    C <- as.numeric(phe.m[phe.m[,"gtsM"] == "C" & as.numeric(phe.m[,1]) >= 20, 1])
    B <- as.numeric(phe.m[phe.m[,"gtsM"] == "B" & as.numeric(phe.m[,1]) >= 20, 1])

    H <- as.numeric(phe.p[phe.p[,"gtsP"] == "H" & as.numeric(phe.p[,1]) >= 20, 1])
    D <- as.numeric(phe.p[phe.p[,"gtsP"] == "D" & as.numeric(phe.p[,1]) >= 20, 1])

    setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny")

      std <- function(x) sd(x)/sqrt(length(x))

      maxD <- max(as.numeric( phe[which(phe[, "sex"] == sex),1]))
      x <- sort(as.numeric( phe[which(phe[, "sex"] == sex),1]))[seq(0,nrow(phe[which(phe[, "sex"] == sex),]), 80)]
      windowM <- cbind(c(0,x), c(x,maxD))
      mids <- c()
      mm <- c()
      for(x in 1:nrow(windowM)){
        day <- windowM[x, 1]
        day2 <- windowM[x, 2]
        mids <- c(mids, (day+day2) / 2)

        nCdR <- (length(which(C > day & C <= day2)) / length(c(C)))
        nBdR <- (length(which(B > day & B <= day2)) / length(c(B)))
        nHdR <- (length(which(H > day & H <= day2)) / length(c(H)))
        nDdR <- (length(which(D > day & D <= day2)) / length(c(D)))
        mm <- rbind(mm, c(nCdR, nBdR, nHdR, nDdR))
      }

      rownames(mm) <- mids
      colnames(mm) <- c("C", "B", "H", "D")
      toR <- which(apply(apply(mm,1, function(x){x == 0}),2, sum) > 0)
      if(length(toR) > 0){
        mm <- mm[-toR,]
        mids <- mids[-toR]
      }
    mp <- mm

    # TODO: Add a test and significance to the plots (ChiSq test)
    # TODO: Standard deviation based on the surrounding bins
    # TODO: Genotype based plots corresponding to the actuarial plots
    # TODO: Paternal colors (blueish) Maternal colors (Redish)
    # TODO: Histogram of all epistasis values with bars of about 0.1
    pdf(paste0("July25/Vita_HazardRatios/", mname, "_",ss[sex],"_Deathrate_FixAxis.pdf"), width = 16.4, height = 12)
      op <- par(cex = 2)
      par(mar = c(5, 5, 4, 3))
      dt1 <- log2(mp[,1]/mp[,2])
      dt2 <- log2(mp[,3]/mp[,4])
      psequence <- mids
      span <- 0.3
      m1 <- loess(dt1 ~ psequence, span = span, degree = 2)
      m2 <- loess(dt2 ~ psequence, span = span, degree = 2)
      smooth1 <- predict(m1)
      smooth2 <- predict(m2)

      E1 <- sum(summary(m1)$residuals^2)
      E2 <- sum(summary(m2)$residuals^2)

      xx <- max(abs(c(dt1, dt2)))

      msex <- "Females"
      if(sex == 1) msex <- "Males"
      plot(c(42, 1350), c(-1.2, 1.2) , t = "n", main = paste0(mname, " ", msex), 
           sub = paste0("Death = 80, alpha: ",span, " E1/E2: ", round(E1,2), "/",round(E2,2)), xaxs="i",
           xaxt="n",yaxt="n", xlab = "Age [d]", ylab = "log2(Hazard Ratio)")

      col1 <- col.main[c(1:2)][1-as.numeric(smooth1 > 0)+1]
      col1[which(smooth1 == 0)] <- "black"

      col1a <- col.main[c(1:2)][1-as.numeric(dt1 > 0)+1]
      col1a[which(dt1 == 0)] <- "black"

      points(psequence, abs(dt1), t = "p", pch=3, col= col1a, xaxt="n")
      points(psequence, abs(smooth1), t = "l", pch=19, col= "black", xaxt="n")
      points(psequence, abs(smooth1), t = "p", pch=19, col= col1, xaxt="n")
      #points(psequence, lods1 * 10, t = "l", xaxt="n")

      col2 <- col.main[c(3:4)][1-as.numeric(smooth2 > 0)+1]
      col2[which(smooth2 == 0)] <- "black"

      col2a <- col.main[c(3:4)][1-as.numeric(dt2 > 0)+1]
      col2a[which(dt2 == 0)] <- "black"

      points(psequence, -abs(dt2), t = "p", pch=3, col= col2a, xaxt="n")
      points(psequence, -abs(smooth2), t = "l", pch=19, col= "black", xaxt="n")
      points(psequence, -abs(smooth2), t = "p", pch=19, col= col2)
      #points(psequence, -(lods2 * 10), t = "l", xaxt="n")

      axis(1, at = seq(0, 1600, 200), seq(0, 1600, 200), las=1)
      axis(1, at = seq(0, 1600, 100), rep("", length(seq(0, 1600, 100))), las=1)
      axis(2, at = seq(-3.6, 3.6, 0.2), round(abs(seq(-3.6, 3.6, 0.2)),1), las=1)
      axis(2, at = seq(-3.6, 3.6, 0.1), rep("", length(abs(seq(-3.6, 3.6, 0.1)))), las=1)
     # axis(2, at = seq(-4, 4, .5), rep("", length(abs(seq(-4, 4, .5)))), las=1)
      text(150, -1.5, "Paternal")
      text(150, 1.5, "Maternal")
      abline(h = 0)
      legend("topright", c("C", "B"), col = col.main[1:2], pch=19, bty = "n")
      legend("bottomright", c("H", "D"), col = col.main[3:4], pch=19, bty = "n")
    dev.off()
  }
}

  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny")
  pdf(paste0(mname, "_RelativeDeaths_15day_WholePopulation.pdf"), width = 16, height = 12)

  par(cex=2)
  par(cex.axis=1.2)


    plot(c(200, 1500), c(0, 1), t = "n", 
         ylab = "cumsum of % deaths in window (relative to population)", xlab = "Lifespan Cut-off Age (days)", yaxt="n", xaxs="i", main = "Effect Curve")
    msequence <- seq(20, 1500, 15)
    std <- function(x) sd(x)/sqrt(length(x))
    mm <- c()
    for(day in msequence){
      day2 <- day + 15
      nCdR <- length(which(C > day & C <= day2)) / length(c(C,D))
      nDdR <- length(which(D > day & D <= day2)) / length(c(C,D))
      mm <- rbind(mm, c(nCdR, nDdR))
    }
    abline(h = seq(0, 15, 1), lty=2, col = "gray")
    points(msequence, cumsum(mm[,1]), t = 'l', col = "#FF3399", lwd=2)
    points(msequence, cumsum(mm[,2]), t = 'l', col = "#CCCC99", lwd=2)

    legend("topleft", c("H", "D"), col = c("#FF3399", "#CCCC99"), pch=18)
    axis(2, at = seq(0, 1, 0.1), seq(0, 1, 0.1), las=2)
    #dev.off()

    abline(v = 970)
    abline(v = 520)
  dev.off()


#  hist(mm[,1] - mm[,2], col = rgb(1,0,0,0.5), breaks = seq(0, 0.02, .0001))
#  hist(mm[,2], add = TRUE, col = rgb(0,1,0,0.5), breaks = seq(0, 0.02, .0001))

  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny")
  pdf(paste0(mname, "_RelativeDeaths_15day_SubPopulation.pdf"), width = 16, height = 12)

  par(cex=2)
  par(cex.axis=1.2)

  dt <- (mm[,1]-mm[,2])

  predict(loess(dt ~ msequence))

  plot(msequence, 1000*smooth(mm[,1]-mm[,2]), t = "b", pch=19, 
       col= c("#FF3399", "#CCCC99")[1-as.numeric(smooth(mm[,1]-mm[,2]) > 0)+1], xaxt="n", ylab = "Relative to mean mortality rate", main = "Vita2b")
  axis(1, at = seq(0, 1600, 200), seq(0, 1600, 200), las=1)
  abline(h = 0)
  legend("topleft", c("H", "D"), col = c("#FF3399", "#CCCC99"), pch=18)

  dev.off()

  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny")
  pdf("Vita2B_Deaths_15day_OwnAllele.pdf", width = 16, height = 12)

  par(cex=2)
  par(cex.axis=1.2)


    plot(c(200, 1500), c(0, 1), t = "n", ylab = "cumsum of % deaths in window (relative to same allele)", xlab = "Lifespan Cut-off Age (days)", yaxt="n", xaxs="i", main = "Effect Curve")
    msequence <- seq(20, 1500, 15)
    std <- function(x) sd(x)/sqrt(length(x))
    mm <- c()
    for(day in msequence){
      day2 <- day + 15
      nCdR <- length(which(C > day & C <= day2)) / length(c(C))
      nDdR <- length(which(D > day & D <= day2)) / length(c(D))
      mm <- rbind(mm, c(nCdR, nDdR))
    }
    abline(h = seq(0, 15, 1), lty=2, col = "gray")
    points(msequence, cumsum(mm[,1]), t = 'l', col = "#FF3399", lwd=2)
    points(msequence, cumsum(mm[,2]), t = 'l', col = "#CCCC99", lwd=2)
    #polygon(c(msequence, rev(msequence)), c(mm[,1] + mm[,2], rev(mm[,1] - mm[,2])), col = rgb(0,0,1,0.5), border = NA)

    #points(msequence, mm[,3], t = 'l', col = rgb(.4,0.2,1,1), lwd=2)
    #polygon(c(msequence, rev(msequence)), c(mm[,3] + mm[,4], rev(mm[,3] - mm[,4])), col = rgb(.4,0.2,1,0.5), border = NA)


    legend("topleft", c("H", "D"), col = c("#FF3399", "#CCCC99"), pch=18)
    axis(2, at = seq(0, 1, 0.1), seq(0, 1, 0.1), las=2)
    #dev.off()

    abline(v = 970)
    abline(v = 520)
  dev.off()

}


