library(RColorBrewer)
library(svglite)
library(vioplot)

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

lodM <- read.table("vita_interactions.txt", sep = "\t")[,-1]

colz.c <- colorRampPalette(brewer.pal(9, "PuRd")[-c(8:9)])(15)

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
pdf(paste0("Figure_3_Interactions.pdf"), width = 24, height = 12)

plot(c(1.5, 30.5), c(3.5, 31.5), t = "n", xaxt='n', yaxt='n', xlab="",ylab="", xaxt="n", yaxt="n",bty="n")
for(x in 1:32){
  for(y in 3:32){
    m1 <- gsub("Vita", "", rownames(lodM)[x])
    m2 <- gsub("Vita", "", colnames(lodM)[y])
    m1 <- substr(m1, 1, nchar(m1)-1)
    m2 <- substr(m2, 1, nchar(m2)-1)
    if(m1 != m2 && !is.na(lodM[x,y])){
      rect(x-0.5,y-0.5, x+0.5, y + 0.5, col = colz.c[round(lodM[x,y])], border = "gray")
      text(x, y, paste0(formatC(lodM[x,y], digits = 1, format = "f"), ""))
    }
  }
}

abline(h = c(3, 6,8,11,12,14,15,19,20,22,23,25,26,27,30,31) - 0.5, lwd=1, lty=2)
abline(v = c(3, 6,8,11,12,14,15,19,20,22,23,25,26,27,30,31) + 0.5, lwd=1, lty=2)

for(x in 1:32){
  for(y in 3:32){
    m1 <- gsub("Vita", "", rownames(lodM)[x])
    m2 <- gsub("Vita", "", colnames(lodM)[y])
    m1 <- substr(m1, 1, nchar(m1)-1)
    m2 <- substr(m2, 1, nchar(m2)-1)
    if(m1 != m2 && !is.na(lodM[x,y])){
      if(lodM[x,y] >= 5) rect(x-0.5,y-0.5, x+0.5, y + 0.5, border = "red",lwd=2)
    }
  }
}


axis(1, at = 1:32, rownames(lodM)[-33],las=2)
axis(2, at = 1:32, rownames(lodM)[-1],las=2)
dev.off()


### Interaction plot
setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

marker1 <- "9_104091597" # 
marker1 <- "6_54992703"
marker2 <- "9_124029281" # 
marker2 <- "11_82178599"

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
col.main <- c("#00A654", "#004BAD", "#B16BE6", "#F02D27")
add.alpha <- function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
col.alpha <- add.alpha(col.main, 0.1)

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
pdf(paste0("Figure_3_Vita6a_vs_Vita11b.pdf"), width = 8, height = 8)

colz <- colorRampPalette(c("#fd8d3c", "#41ab5d"))(10)
plot(c(0.5,4.5), c(0.5, 4.5), t = "n", main = "Interaction Vita6a & Vita11b (>=365 days)", 
     xlab = "Haplotype Vita6a", ylab = "Haplotype Vita11b", xaxt="n", yaxt = "n", yaxs= "i", xaxs= "i")
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
plot(c(1,4), c(820, 900), t = "n", main = "Interaction Vita6a & Vita11b (>=365 days)", 
     xlab = "Haplotype Vita6a", ylab = "Adjusted Longevity (days)", xaxt="n")
off <- -0.15
j <- 1
for(v1 in c("AC", "AD", "BC", "BD")){
  i <- 1
  vals <- c()
  errs <- c()
  for(v2 in c("AC", "AD", "BC", "BD")){
    values <- cdata[which(cdata[, "gts1"] == v1 & cdata[, "gts2"] == v2), "pheAdj"]
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
legend("bottomleft", c("CH", "CD", "BH", "BD"), col = col.main, lwd=1, pch = 20, title = "Vita11b")

dev.off()



