#
# Figure3.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Figures for the main paper text (Main text Figure 3, but some figure might have moved into other figures to make them fit)
# Figures: 
# - G x G Interaction split by males and females
# - G x G Interaction figure of Vita4a vs Vita10a
# - Epistasis by Sex visualization showing the difference in epistasis between males and females (now in Figure 5)
# - G x G Interaction figures of all Vita loci versus all Vita loci (Figure 5 and Supplements)
# - G x G Interaction figures of the green dots seen in figure 5 (shared between males and females but not significant in either
#

library(RColorBrewer)
library(svglite)
library(vioplot)

all <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085")

names(all) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa")

lodM.m <- read.table(paste0("DataSet/output/vita_soma_interactions_2way_males_tp42.txt"), sep = "\t")
lodM.f <- read.table(paste0("DataSet/output/vita_soma_interactions_2way_females_tp42.txt"), sep = "\t")

colz.c <- colorRampPalette(c("white", "lightskyblue3"))(15)
colz.c2 <- colorRampPalette(c("white", "plum2"))(15)

pdf(paste0("DataSet/output/Figure_4_GxG_Interaction_combined.pdf"), width = 24, height = 12)
plot(c(1.5, 29.5), c(1.5, 29.5), t = "n", xaxt='n', yaxt='n', xlab="",ylab="", xaxt="n", yaxt="n",bty="n")
axis(1, at = 1:29, rownames(lodM.m),las=2)
axis(2, at = 1:29, rev(rownames(lodM.m)),las=2)
for(x in 1:30){
  for(y in 1:30){
    xp <- x
    yp <- 31 - y
    m1 <- gsub("Vita", "", rownames(lodM.m)[x])
    m2 <- gsub("Vita", "", colnames(lodM.m)[y])
    m1 <- substr(m1, 1, nchar(m1)-1)
    m2 <- substr(m2, 1, nchar(m2)-1)
    if(x > y && !is.na(lodM.m[x,y])){
      rect(xp-0.5,yp-0.5, xp+0.5, yp + 0.5, col = colz.c[1+round(lodM.m[x,y])], border = "gray")
      text(xp, yp, paste0(formatC(lodM.m[x,y], digits = 1, format = "f"), ""))
    }
    if(x < y &&  !is.na(lodM.f[y,x])){
      rect(xp-0.5,yp-0.5, xp+0.5, yp + 0.5, col = colz.c2[1+round(lodM.f[y,x])], border = "gray")
      text(xp, yp, paste0(formatC(lodM.f[y,x], digits = 1, format = "f"), ""))
    }
  }
}
box()

#abline(h = c(4,7,9,10,11,13,16,17,19,20,21,22,23,24,25,27) - 0.5, lwd=1, lty=2)
#abline(v = c(3,6,8,9,10,12,15,16,18,19,20,21,22,23,24,26) + 0.5, lwd=1, lty=2)

for(x in 1:29){
  for(y in 1:29){
    xp <- x
    yp <- 31 - y
    m1 <- gsub("Vita", "", rownames(lodM.m)[x])
    m2 <- gsub("Vita", "", colnames(lodM.m)[y])
    m1 <- substr(m1, 1, nchar(m1)-1)
    m2 <- substr(m2, 1, nchar(m2)-1)
    if(!is.na(lodM.m[x,y])){
    #  if(x <= y && lodM.m[x,y] >= 4) rect(xp-0.5,yp-0.5, xp+0.5, yp + 0.5, border = "red",lwd=2)
    }
    if(!is.na(lodM.f[y,x])){
     # if(x >= y && lodM.f[y,x] >= 4) rect(xp-0.5,yp-0.5, xp+0.5, yp + 0.5, border = "red",lwd=2)
    }
  }
}

#dev.off()


### Interaction plot
library(qtl)

source("ActuarialMapping/adjustXprobs.R")

# Read cross object
mcross <- read.cross(format="csvr", file="DataSet/um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

marker1 <- "4_52524395" # 
marker2 <- "10_72780332" # 

mp1 <- gtsp[, grep(marker1, colnames(gtsp))]
gts1 <- unlist(lapply(lapply(lapply(apply(mp1,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
  if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
}))

mp2 <- gtsp[, grep(marker2, colnames(gtsp))]
gts2 <- unlist(lapply(lapply(lapply(apply(mp2,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
  if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
}))


cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]),
                    gts1 = gts1, gts2 = gts2)

idx <- which(cdata[, "longevity"] >= 365)
cdata <- cdata[idx,]

mlm <- lm("longevity ~ sex + site + cohort + treatment", data = cdata)

pheAdj <- rep(NA, nrow(cdata))
names(pheAdj) <-  rownames(cdata)

adj <- residuals(mlm) + mean(cdata[, "longevity"])
pheAdj[names(adj)] <- round(adj)

cdata <- cbind(pheAdj, cdata)
cdata <- cdata[which(!is.na(cdata[, "gts1"]) & !is.na(cdata[, "gts2"])),]
std <- function(x) sd(x)/sqrt(length(x))
col.main <- c("#01A654", "#1750A3", "#714F99", "#EE3129")
add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
col.alpha <- add.alpha(col.main, 0.1)
col.alpha2 <- add.alpha(col.main, 0.3)

pdf(paste0("DataSet/output/Figure_3_Vita4a_vs_Vita10a.pdf"), width = 8, height = 8)

colz <- colorRampPalette(c("#fd8d3c", "#41ab5d"))(10)
plot(c(0.5,4.5), c(0.5, 4.5), t = "n", main = "Interaction Vita4a & Vita10a (>=365 days)", 
     xlab = "Haplotype Vita4a", ylab = "Haplotype Vita10a", xaxt="n", yaxt = "n", yaxs= "i", xaxs= "i")
x <- 1
for(v1 in c("AC", "AD", "BC", "BD")){
  y <- 1
  for(v2 in c("AC", "AD", "BC", "BD")){
    values <- cdata[which(cdata[, "gts1"] == v1 & cdata[, "gts2"] == v2), "pheAdj"]
    cat(log(mean(values)), "\n")
    cI <- (mean(values) - 800) %/% 10
    points(x, y, cex = cI * 3.0, pch=20, col = colz[cI])
    text(x, y, round(mean(values)), col = "black")
    text(x, y-0.1, paste0("n = ", round(length(values))), col = "black", cex=0.7)
    y <- y+1
  }
  x <- x+1
}
abline(h=seq(0.5, 4.5, 1))
abline(v=seq(0.5, 4.5, 1))
axis(2, at = 1:4, c("CH", "CD", "BH", "BD"), las=2)
axis(1, at = 1:4, c("CH", "CD", "BH", "BD"))
dev.off()

### OLD

### TODO: Vita1b & Vita9a
### TODO: Vita1c & Vita3b

all <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085")

names(all) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa")

lodM.m <- read.table(paste0("DataSet/output/vita_interactions_2way_males_tp42.txt"), sep = "\t")
lodM.f <- read.table(paste0("DataSet/output/vita_interactions_2way_females_tp42.txt"), sep = "\t")

males <- unlist(lodM.m[lower.tri(lodM.m)])
females <-unlist(lodM.f[lower.tri(lodM.f)])

same <- (abs(males - females) < 0.5)
m <- (males - females > 0.5)
f <- (females - males > 0.5)

colz <- rep("#00A651", length(males))
colz[which(m)] <- "#00AEEF"
colz[which(f)] <- "#FF3333"
colz[which((males^2 + females^2) < 5)] <- "black"

ii <- which(colz == "#00A651" | colz == "#00AEEF" | colz == "#FF3333"| colz == "black")
idx <- sort(males[ii], dec = TRUE, index.return = TRUE)$ix
mmm <- males[ii][idx]

pdf(paste0("DataSet/output/SexByEpistasis_Cor_TP42_April25.pdf"), width = 24, height = 24)

plot(c(0,6), c(0,6), pch = 19, col = colz, t = "n", 
     xlab = "LOD females", ylab = "LOD males", main = "GxG correlation between LODs")

abline(h = 3.8, lty=2)
abline(v = 3.8, lty=2)
i <- 1
for(x in mmm){
  p1 <- names(which(apply(lodM.m == x,1,any, na.rm=TRUE)))
  p2 <- names(which(apply(lodM.m == x,2,any, na.rm=TRUE)))
  chr1 <- substring(gsub("Vita", "", p1), 1, nchar(gsub("Vita", "", p1))-1)
  chr2 <- substring(gsub("Vita", "", p2), 1, nchar(gsub("Vita", "", p2))-1)
  if(chr1 != chr2) {
    col <- "black"
    if(abs(lodM.f[p1,p2] - lodM.m[p1,p2]) < 0.5) col <- "#00A651"
    if(lodM.f[p1,p2] - lodM.m[p1,p2] > 0.5) col <- "#FF3333"
    if(lodM.m[p1,p2] - lodM.f[p1,p2] > 0.5) col <- "#00AEEF"
    if((lodM.m[p1,p2]^2 + lodM.f[p1,p2]^2) < 5) col <- "black"
    points(lodM.f[p1,p2], lodM.m[p1,p2], col = col, pch = 19)
    if(col != "black") text(lodM.f[p1,p2], lodM.m[p1,p2], gsub("Vita", "", paste0(p1, "-",p2)))
  }
  i <- i + 1
  #cat(p1, " ", p2, "=", x,"\n")
}
dev.off()

### GxG plots for all SIGNIFICANT GxG interaction pairs

cor.test(unlist(lodM.m[lower.tri(lodM.m)])[ii], unlist(lodM.f[lower.tri(lodM.f)])[ii], method = "spearman")

for(timepoint in c(42, 365, 740, 905)){
  for(mX1 in 1:(length(all)-1)){
    for(mX2 in mX1:length(all)){
      m1 <- names(all)[mX1]
      m2 <- names(all)[mX2]
      if(m1 == m2) next;

      if(lodM.m[m2,m1] >= 3.8 || lodM.f[m2,m1] >= 3.8){

      pdf(paste0("DataSet/output/GxG_",timepoint,"_",m1,"_",m2, "_new.pdf"), width = 18, height = 8)

      marker1 <- all[m1]
      marker2 <- all[m2]


      mp1 <- gtsp[, grep(marker1, colnames(gtsp))]
      gts1 <- unlist(lapply(lapply(lapply(apply(mp1,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
        if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
      }))

      mp2 <- gtsp[, grep(marker2, colnames(gtsp))]
      gts2 <- unlist(lapply(lapply(lapply(apply(mp2,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
        if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
      }))

      cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                          sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                          site = as.factor(pull.pheno(mcross)[, "site"]),
                          cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                          treatment = as.factor(pull.pheno(mcross)[, "treatment"]),
                          gts1 = gts1, gts2 = gts2)

      cdata <- cdata[which(cdata[, "longevity"] > timepoint),]
      cdata <- cdata[which(!is.na(cdata[, "gts1"]) & !is.na(cdata[, "gts2"])),]

      mSm <- rbind(c(0, 20, 1),
        c(21, 50, 2),
        c(51, 100, 3),
        c(101, 150, 4),
        c(151, 500, 5))

      pSize <- function(x){ return(mSm[mSm[,1] < x & mSm[,2] > x,3]) }

      op <- par(mfrow = c(1,3))
      par(cex=4)
      par(cex.axis=2)
      par(cex.main=2)
      par(cex.sub=2)
      par(cex.lab=2)


      ### Combined
      op <- par(mfrow = c(1,3))
      op <- par(mar = c(5.1, 5.1, 4.1, 0.1))
      xdata <- cdata

      lodMS <- read.table(paste0("DataSet/output/vita_interactions_2way_combined_tp",timepoint,".txt"), sep = "\t")

      plot(c(0.5, 4.5), c(-80, +80), t = "n", 
           main = paste0("Interaction ",m1," & ",m2," (Combined)", paste0("\nT",timepoint,", LOD = ",round(lodMS[m2,m1],1))),
           xlab = paste0(m1, " Genotype"), ylab = "Expectancy (d)", xaxt="n", yaxt="n", yaxs = "i")
      text(x=1.4, y = 71, labels = bquote(atop(bold(bolditalic(.(m1)) ~ x ~ bolditalic(.(m2))), bold(at ~ T[.(timepoint)]))), cex=2)
      text(x=3.75, y = 75, labels = bquote(Combined), cex=2)
      rect(0.75, -80, 1.25, 70, col = col.alpha2[1], border = NA)
      rect(1.75, -80, 2.25, 70, col = col.alpha2[2], border = NA)
      rect(2.75, -80, 3.25, 70, col = col.alpha2[3], border = NA)
      rect(3.75, -80, 4.25, 70, col = col.alpha2[4], border = NA)

      off <- -0.15
      j <- 1
      minNs <- c()
      for(v1 in c("AC", "AD", "BC", "BD")){
        i <- 1
        vals <- c()
        errs <- c()
        sizes <- c()
        v1M <- xdata[which(xdata[, "gts1"] == v1), "longevity"]
        for(v2 in c("AC", "AD", "BC", "BD")){
          values <- xdata[which(xdata[, "gts1"] == v1 & xdata[, "gts2"] == v2), "longevity"]
          #vioplot(values, at = i + off, add = TRUE, wex=0.1, col = col.main[j])
          vals <- c(vals, mean(values))
          errs <- c(errs, std(values))
          sizes <- c(sizes, pSize(length(values)))
          minNs <- c(minNs, length(values))
          cat(m1, " ", m2, "size:", length(values), "\n")
          i <- i + 1
        }
        min(minNs,na.rm=TRUE)
        points(off + 1:4, vals - mean(v1M), col=col.main[j], t = "b",pch=20, cex=sizes)

        points(off + 1:4, vals - mean(v1M) + errs, col=col.main[j], t = "p",pch="-")
        points(off + 1:4, vals - mean(v1M) - errs, col=col.main[j], t = "p",pch="-")
        for(x in 1:4){
          points(c(off + (1:4)[x], off + (1:4)[x]), c(vals[x] - mean(v1M) - errs[x], vals[x] - mean(v1M) + errs[x]), col=col.main[j], t = "l")
        }

        #polygon(c(1:4, 4:1), c(vals - mean(v1M) + errs, rev(vals - mean(v1M) - errs)), col = col.alpha[j], border = NA)
        off <- off + 0.1
        j <- j + 1
      }
      axis(1, at = 1, bquote(bold("CH")), col.axis= col.main[1])
      axis(1, at = 2, bquote(bold("CD")), col.axis= col.main[2])
      axis(1, at = 3, bquote(bold("BH")), col.axis= col.main[3])
      axis(1, at = 4, bquote(bold("BD")), col.axis= col.main[4])
      axis(2, at = seq(-80,80,20), seq(-80,80,20), las=2)
      legend("bottomleft", c("CH", "CD", "BH", "BD"), col = col.main, lwd=1, pch = 20, title = m2, cex=1.5,bg = "white")
      legend("bottomright", c("0-20", "21-50", "51-100", "101-150", ">=151"), pch = 20, pt.cex = 1:5, title = paste0("minN=", min(minNs,na.rm=TRUE)), cex=1.5,bg = "white")

      ### Female
      xdata <- cdata[which(cdata[, "sex"] == 0),]

      lodMS <- read.table(paste0("DataSet/output/vita_interactions_2way_females_tp",timepoint,".txt"), sep = "\t")
      plot(c(0.5, 4.5), c(-80, 80), t = "n", 
           main = paste0("Interaction ",m1," & ",m2," (Females)", paste0("\nT",timepoint,", LOD = ",round(lodMS[m2,m1],1))),
           xlab = paste0(m1, " Genotype"), ylab = "", xaxt="n", yaxt="n", yaxs = "i")
      text(x=1.4, y = 71, labels = bquote(atop(bold(bolditalic(.(m1)) ~ x ~ bolditalic(.(m2))), bold(at ~ T[.(timepoint)]))), cex=2)
      text(x=3.8, y = 75, labels = bquote(Females), cex=2)
      rect(0.75, -80, 1.25, 70, col = col.alpha2[1], border = NA)
      rect(1.75, -80, 2.25, 70, col = col.alpha2[2], border = NA)
      rect(2.75, -80, 3.25, 70, col = col.alpha2[3], border = NA)
      rect(3.75, -80, 4.25, 70, col = col.alpha2[4], border = NA)

      off <- -0.15
      j <- 1
      minNs <- c()
      for(v1 in c("AC", "AD", "BC", "BD")){
        i <- 1
        vals <- c()
        errs <- c()
        sizes <- c()
        v1M <- xdata[which(xdata[, "gts1"] == v1), "longevity"]
        for(v2 in c("AC", "AD", "BC", "BD")){
          values <- xdata[which(xdata[, "gts1"] == v1 & xdata[, "gts2"] == v2), "longevity"]
          #vioplot(values, at = i + off, add = TRUE, wex=0.1, col = col.main[j])
          vals <- c(vals, mean(values))
          errs <- c(errs, std(values))
          sizes <- c(sizes, pSize(length(values)))
          minNs <- c(minNs, length(values))
          cat(m1, " ", m2, "size:", length(values), "\n")
          i <- i + 1
        }
        points(off + 1:4, vals - mean(v1M), col=col.main[j], t = "b",pch=20, cex=sizes)

        points(off + 1:4, vals - mean(v1M) + errs, col=col.main[j], t = "p",pch="-")
        points(off + 1:4, vals - mean(v1M) - errs, col=col.main[j], t = "p",pch="-")
        for(x in 1:4){
          points(c(off + (1:4)[x], off + (1:4)[x]), c(vals[x] - mean(v1M) - errs[x], vals[x] - mean(v1M) + errs[x]), col=col.main[j], t = "l")
        }

        #polygon(c(1:4, 4:1), c(vals - mean(v1M) + errs, rev(vals - mean(v1M) - errs)), col = col.alpha[j], border = NA)
        off <- off + 0.1
        j <- j + 1
      }


      axis(1, at = 1, bquote(bold("CH")), col.axis= col.main[1])
      axis(1, at = 2, bquote(bold("CD")), col.axis= col.main[2])
      axis(1, at = 3, bquote(bold("BH")), col.axis= col.main[3])
      axis(1, at = 4, bquote(bold("BD")), col.axis= col.main[4])
      legend("bottomright", c("0-20", "21-50", "51-100", "101-150", ">=151"), pch = 20, pt.cex = 1:5, title = paste0("minN=", min(minNs,na.rm=TRUE)), cex=1.5,bg = "white")

      axis(2, at = seq(-80,80,20), rep("",9), las=2)
      #legend("bottomleft", c("CH", "CD", "BH", "BD"), col = col.main, lwd=1, pch = 20, title = m2)

      ## Males
      xdata <- cdata[which(cdata[, "sex"] == 1),]

      lodMS <- read.table(paste0("DataSet/output/vita_interactions_2way_males_tp",timepoint,".txt"), sep = "\t")
      plot(c(0.5, 4.5), c(-80, 80), t = "n", 
           main = paste0("Interaction ",m1," & ",m2," (Males)", paste0("\nT",timepoint,", LOD = ",round(lodMS[m2,m1],1))),
           xlab = paste0(m1, " Genotype"), ylab = "", xaxt="n", yaxt="n", yaxs = "i")
      text(x=1.4, y = 71, labels = bquote(atop(bold(bolditalic(.(m1)) ~ x ~ bolditalic(.(m2))), bold(at ~ T[.(timepoint)]))), cex=2)
      text(x=3.85, y = 75, labels = bquote(Males), cex=2)
      rect(0.75, -80, 1.25, 70, col = col.alpha2[1], border = NA)
      rect(1.75, -80, 2.25, 70, col = col.alpha2[2], border = NA)
      rect(2.75, -80, 3.25, 70, col = col.alpha2[3], border = NA)
      rect(3.75, -80, 4.25, 70, col = col.alpha2[4], border = NA)

      off <- -0.15
      j <- 1
      minNs <- c()
      for(v1 in c("AC", "AD", "BC", "BD")){
        i <- 1
        vals <- c()
        errs <- c()
        sizes <- c()
        v1M <- xdata[which(xdata[, "gts1"] == v1), "longevity"]
        for(v2 in c("AC", "AD", "BC", "BD")){
          values <- xdata[which(xdata[, "gts1"] == v1 & xdata[, "gts2"] == v2), "longevity"]
          #vioplot(values, at = i + off, add = TRUE, wex=0.1, col = col.main[j])
          vals <- c(vals, mean(values))
          errs <- c(errs, std(values))
          sizes <- c(sizes, pSize(length(values)))
          minNs <- c(minNs, length(values))
          i <- i + 1
        }
        points(off + 1:4, vals - mean(v1M), col=col.main[j], t = "b",pch=20, cex=sizes)

        points(off + 1:4, vals - mean(v1M) + errs, col=col.main[j], t = "p",pch="-")
        points(off + 1:4, vals - mean(v1M) - errs, col=col.main[j], t = "p",pch="-")
        for(x in 1:4){
          points(c(off + (1:4)[x], off + (1:4)[x]), c(vals[x] - mean(v1M) - errs[x], vals[x] - mean(v1M) + errs[x]), col=col.main[j], t = "l")
        }

        #polygon(c(1:4, 4:1), c(vals - mean(v1M) + errs, rev(vals - mean(v1M) - errs)), col = col.alpha[j], border = NA)
        off <- off + 0.1
        j <- j + 1
      }
      axis(1, at = 1, bquote(bold("CH")), col.axis= col.main[1])
      axis(1, at = 2, bquote(bold("CD")), col.axis= col.main[2])
      axis(1, at = 3, bquote(bold("BH")), col.axis= col.main[3])
      axis(1, at = 4, bquote(bold("BD")), col.axis= col.main[4])

      legend("bottomright", c("0-20", "21-50", "51-100", "101-150", ">=151"), pch = 20, pt.cex = 1:5, title = paste0("minN=", min(minNs,na.rm=TRUE)), cex=1.5,bg = "white")

      axis(2, at = seq(-80,80,20), rep("",9), las=2)

      dev.off()
      }
    }
  }
}

### DO the green dots
combos <- rbind(
            c("Vita2a", "Vita17a"),
            c("Vita11a", "Vita14a"),
            c("Vita4a", "Vita10a"),
            c("Vita10a", "Vita11a"),
            c("Vita3a", "Vita6b"),
            c("Vita11b", "Vita13a"),
            c("Vita11a", "Vita13a"),
            c("Vita11a", "Vita15a"),
            c("Vita2c", "Vita18a"),
            c("Vita3b", "Vita4a"),
            c("Vita1b", "Vita10a"),
            c("Vita1c", "Vita13a"),
            c("Vita6a", "Vita12a"))


for(timepoint in c(42, 365, 740, 905)){
  for(pp in 1:nrow(combos)){
      m1 = combos[pp,1]
      m2 = combos[pp,2]

      pdf(paste0("DataSet/output/Green_Dot_GxG_",timepoint,"_",m1,"_",m2, ".pdf"), width = 18, height = 8)

      marker1 <- all[m1]
      marker2 <- all[m2]


      mp1 <- gtsp[, grep(marker1, colnames(gtsp))]
      gts1 <- unlist(lapply(lapply(lapply(apply(mp1,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
        if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
      }))

      mp2 <- gtsp[, grep(marker2, colnames(gtsp))]
      gts2 <- unlist(lapply(lapply(lapply(apply(mp2,1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
        if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
      }))

      cdata <- data.frame(longevity = pull.pheno(mcross)[, "longevity"], 
                          sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                          site = as.factor(pull.pheno(mcross)[, "site"]),
                          cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                          treatment = as.factor(pull.pheno(mcross)[, "treatment"]),
                          gts1 = gts1, gts2 = gts2)

      cdata <- cdata[which(cdata[, "longevity"] > timepoint),]
      cdata <- cdata[which(!is.na(cdata[, "gts1"]) & !is.na(cdata[, "gts2"])),]

      mSm <- rbind(c(0, 20, 1),
        c(21, 50, 2),
        c(51, 100, 3),
        c(101, 150, 4),
        c(151, 500, 5))

      pSize <- function(x){ return(mSm[mSm[,1] < x & mSm[,2] > x,3]) }

      op <- par(mfrow = c(1,3))
      par(cex=4)
      par(cex.axis=2)
      par(cex.main=2)
      par(cex.sub=2)
      par(cex.lab=1.5)
      ### Combined
      op <- par(mfrow = c(1,3))
      op <- par(mar = c(5.1, 4.1, 4.1, 0.1))
      xdata <- cdata

      lodMS <- read.table(paste0("DataSet/output/vita_interactions_2way_combined_tp",timepoint,".txt"), sep = "\t")

      plot(c(0.5, 4.5), c(-80, +80), t = "n", 
           main = paste0("Interaction ",m1," & ",m2," (Combined)", paste0("\nT",timepoint,", LOD = ",round(lodMS[m2,m1],1))),
           xlab = paste0("Haplotype ",m1), ylab = "", xaxt="n", yaxt="n", yaxs = "i")
      text(x=1.4, y = 71, labels = bquote(atop(bold(bolditalic(.(m1)) ~ x ~ bolditalic(.(m2))), bold(at ~ T[.(timepoint)]))), cex=2)
      text(x=3.75, y = 75, labels = bquote(Combined), cex=2)
      rect(0.95, -80, 1.05, 70, col = col.alpha2[1], border = NA)
      rect(1.95, -80, 2.05, 70, col = col.alpha2[2], border = NA)
      rect(2.95, -80, 3.05, 70, col = col.alpha2[3], border = NA)
      rect(3.95, -80, 4.05, 70, col = col.alpha2[4], border = NA)

      off <- -0.15
      j <- 1
      for(v1 in c("AC", "AD", "BC", "BD")){
        i <- 1
        vals <- c()
        errs <- c()
        sizes <- c()
        v1M <- xdata[which(xdata[, "gts1"] == v1), "longevity"]
        for(v2 in c("AC", "AD", "BC", "BD")){
          values <- xdata[which(xdata[, "gts1"] == v1 & xdata[, "gts2"] == v2), "longevity"]
          #vioplot(values, at = i + off, add = TRUE, wex=0.1, col = col.main[j])
          vals <- c(vals, mean(values))
          errs <- c(errs, std(values))
          sizes <- c(sizes, pSize(length(values)))
          cat(m1, " ", m2, "size:", length(values), "\n")
          i <- i + 1
        }
        points(vals - mean(v1M)+c(-0.1,-0.05,0.05,0.1), col=col.main[j], t = "b",pch=20, cex=sizes)
        polygon(c(1:4, 4:1), c(vals - mean(v1M) + errs, rev(vals - mean(v1M) - errs)), col = col.alpha[j], border = NA)
        off <- off + 0.1
        j <- j + 1
      }
      axis(1, at = 1, bquote(bold("CH")), col.axis= col.main[1])
      axis(1, at = 2, bquote(bold("CD")), col.axis= col.main[2])
      axis(1, at = 3, bquote(bold("BH")), col.axis= col.main[3])
      axis(1, at = 4, bquote(bold("BD")), col.axis= col.main[4])
      axis(2, at = seq(-80,80,20), seq(-80,80,20), las=2)
      legend("bottomleft", c("CH", "CD", "BH", "BD"), col = col.main, lwd=1, pch = 20, title = m2, cex=1.5,bg = "white")
      legend("bottomright", c("0-20", "21-50", "51-100", "101-150", ">=151"), pch = 20, pt.cex = 1:5, cex=1.5,bg = "white")

      ### Female
      xdata <- cdata[which(cdata[, "sex"] == 0),]
      lodMS <- read.table(paste0("DataSet/output/vita_interactions_2way_females_tp",timepoint,".txt"), sep = "\t")
      plot(c(0.5, 4.5), c(-80, 80), t = "n", 
           main = paste0("Interaction ",m1," & ",m2," (Females)", paste0("\nT",timepoint,", LOD = ",round(lodMS[m2,m1],1))),
           xlab = paste0("Haplotype ",m1), ylab = "", xaxt="n", yaxt="n", yaxs = "i")
      text(x=1.4, y = 71, labels = bquote(atop(bold(bolditalic(.(m1)) ~ x ~ bolditalic(.(m2))), bold(at ~ T[.(timepoint)]))), cex=2)
      text(x=3.8, y = 75, labels = bquote(Females), cex=2)
      rect(0.95, -80, 1.05, 70, col = col.alpha2[1], border = NA)
      rect(1.95, -80, 2.05, 70, col = col.alpha2[2], border = NA)
      rect(2.95, -80, 3.05, 70, col = col.alpha2[3], border = NA)
      rect(3.95, -80, 4.05, 70, col = col.alpha2[4], border = NA)

      off <- -0.15
      j <- 1
      for(v1 in c("AC", "AD", "BC", "BD")){
        i <- 1
        vals <- c()
        errs <- c()
        sizes <- c()
        v1M <- xdata[which(xdata[, "gts1"] == v1), "longevity"]
        for(v2 in c("AC", "AD", "BC", "BD")){
          values <- xdata[which(xdata[, "gts1"] == v1 & xdata[, "gts2"] == v2), "longevity"]
          #vioplot(values, at = i + off, add = TRUE, wex=0.1, col = col.main[j])
          vals <- c(vals, mean(values))
          errs <- c(errs, std(values))
          sizes <- c(sizes, pSize(length(values)))
          cat(m1, " ", m2, "size:", length(values), "\n")
          i <- i + 1
        }
        points(vals - mean(v1M), col=col.main[j], t = "b",pch=20, cex=sizes)
        polygon(c(1:4, 4:1), c(vals - mean(v1M) + errs, rev(vals - mean(v1M) - errs)), col = col.alpha[j], border = NA)
        off <- off + 0.1
        j <- j + 1
      }


      axis(1, at = 1, bquote(bold("CH")), col.axis= col.main[1])
      axis(1, at = 2, bquote(bold("CD")), col.axis= col.main[2])
      axis(1, at = 3, bquote(bold("BH")), col.axis= col.main[3])
      axis(1, at = 4, bquote(bold("BD")), col.axis= col.main[4])

      axis(2, at = seq(-80,80,20), rep("",9), las=2)
      #legend("bottomleft", c("CH", "CD", "BH", "BD"), col = col.main, lwd=1, pch = 20, title = m2)

      ## Males
      xdata <- cdata[which(cdata[, "sex"] == 1),]

      lodMS <- read.table(paste0("DataSet/output/vita_interactions_2way_males_tp",timepoint,".txt"), sep = "\t")
      plot(c(0.5, 4.5), c(-80, 80), t = "n", 
           main = paste0("Interaction ",m1," & ",m2," (Males)", paste0("\nT",timepoint,", LOD = ",round(lodMS[m2,m1],1))),
           xlab = paste0("Haplotype ",m1), ylab = "", xaxt="n", yaxt="n", yaxs = "i")
      text(x=1.4, y = 71, labels = bquote(atop(bold(bolditalic(.(m1)) ~ x ~ bolditalic(.(m2))), bold(at ~ T[.(timepoint)]))), cex=2)
      text(x=3.85, y = 75, labels = bquote(Males), cex=2)
      rect(0.95, -80, 1.05, 70, col = col.alpha2[1], border = NA)
      rect(1.95, -80, 2.05, 70, col = col.alpha2[2], border = NA)
      rect(2.95, -80, 3.05, 70, col = col.alpha2[3], border = NA)
      rect(3.95, -80, 4.05, 70, col = col.alpha2[4], border = NA)

      off <- -0.15
      j <- 1
      for(v1 in c("AC", "AD", "BC", "BD")){
        i <- 1
        vals <- c()
        errs <- c()
        sizes <- c()
        v1M <- xdata[which(xdata[, "gts1"] == v1), "longevity"]
        for(v2 in c("AC", "AD", "BC", "BD")){
          values <- xdata[which(xdata[, "gts1"] == v1 & xdata[, "gts2"] == v2), "longevity"]
          #vioplot(values, at = i + off, add = TRUE, wex=0.1, col = col.main[j])
          vals <- c(vals, mean(values))
          errs <- c(errs, std(values))
          sizes <- c(sizes, pSize(length(values)))
          i <- i + 1
        }
        points(vals - mean(v1M), col=col.main[j], t = "b",pch=20, cex=sizes)
        polygon(c(1:4, 4:1), c(vals - mean(v1M) + errs, rev(vals - mean(v1M) - errs)), col = col.alpha[j], border = NA)
        off <- off + 0.1
        j <- j + 1
      }
      axis(1, at = 1, bquote(bold("CH")), col.axis= col.main[1])
      axis(1, at = 2, bquote(bold("CD")), col.axis= col.main[2])
      axis(1, at = 3, bquote(bold("BH")), col.axis= col.main[3])
      axis(1, at = 4, bquote(bold("BD")), col.axis= col.main[4])

      axis(2, at = seq(-80,80,20), rep("",9), las=2)

      dev.off()
  }
}

