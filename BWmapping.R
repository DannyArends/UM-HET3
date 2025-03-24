setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

#Sample names
snames <- as.character(pull.pheno(mcross)[, "GenoID"])

bwA <- read.csv("bw_RichM.txt", sep = "\t", header = TRUE, comment.char = "#", row.names=1, na.strings = c("NA", "", "x"))

#42 days
bw <- read.csv("ITP_50601.csv", header = TRUE, comment.char = "#", skip=11, row.names=2,na.strings = c("NA", "", "x"))
bw <- bw[which(bw[, "DA2024"] == 1),]

#12 months
trait <- read.csv("ITP_10003.csv", header = TRUE, comment.char = "#", skip=10, row.names=2,na.strings = c("NA", "", "x"))
trait <- trait[which(trait[, "DA2024"] == 1),]

gtsp <- pull.genoprob(mcross)
cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    adjLongevity = NA,
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw42 = as.numeric(bw[snames,"Value"]),
                    adjBw42 = NA,
                    bw6 = as.numeric(pull.pheno(mcross)[, "bw6"]),
                    bw12 = as.numeric(trait[snames, "Value"]),
                    bw18 = as.numeric(pull.pheno(mcross)[, "bw18"]),
                    bw24 = as.numeric(pull.pheno(mcross)[, "bw24"]),
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))
rownames(cdata) <- snames

noData <- snames[which(is.na(cdata[, "bw42"]))]
addData <- noData[which(noData %in% rownames(bwA))]

#Add the data from Rich to the data we had
cdata[addData, "bw42"] <- bwA[addData, "W42d"]

# Our Progressive Mapping Sequence
markers <- unique(unlist(lapply(strsplit(colnames(gtsp), ":"), "[",1)))

lods.f <- c()
for(x in c("bw42", "bw6", "bw12", "bw18", "bw24")){
  idx <- which(cdata[, "sex"] == 0)
  adata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  lods.c <- c()
  lm.null <- lm(adata[,x] ~ site + cohort + treatment + 0, data = adata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(adata[,x] ~ site + cohort + treatment + mp + 0, data = adata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.f  <- rbind(lods.f , lods.c)
  cat("Done", x, "\n")
}
colnames(lods.f) <- colnames(pull.geno(mcross))
rownames(lods.f) <- c("bw42", "bw6", "bw12", "bw18", "bw24")


lods.m <- c()
for(x in c("bw42", "bw6", "bw12", "bw18", "bw24")){
  idx <- which(cdata[, "sex"] == 1)
  adata <- cdata[idx,]
  gtsM <- gtsp[idx,]

  lods.c <- c()
  lm.null <- lm(adata[,x] ~ site + cohort + treatment + 0, data = adata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(adata[,x] ~ site + cohort + treatment + mp + 0, data = adata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.m  <- rbind(lods.m, lods.c)
  cat("Done", x, "\n")
}
colnames(lods.m) <- colnames(pull.geno(mcross))
rownames(lods.m) <- c("bw42", "bw6", "bw12", "bw18", "bw24")


lods.C <- c()
for(x in c("bw42", "bw6", "bw12", "bw18", "bw24")){
  adata <- cdata
  gtsM <- gtsp

  lods.c <- c()
  lm.null <- lm(adata[,x] ~ sex + site + cohort + treatment + 0, data = adata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsM[, grep(marker, colnames(gtsM))]
    lm.alt <- lm(adata[,x] ~ sex + site + cohort + treatment + mp + 0, data = adata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  lods.C  <- rbind(lods.C, lods.c)
  cat("Done", x, "\n")
}
colnames(lods.C) <- colnames(pull.geno(mcross))
rownames(lods.C) <- c("bw42", "bw6", "bw12", "bw18", "bw24")


# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))


subset <- map
subset <- cbind(subset, cpos = NA)
gap <- 10000000
chr.start <- c(0)
chr.mids <- c()
cp <- 0
y <- 1
for(x in c(1:19, "X")){
  cl <- max(as.numeric(subset[which(subset[,1] == x),2]))
  chr.start <- c(chr.start, cl + cp + gap)
  subset[which(subset[,1] == x), "cpos"] <- as.numeric(subset[which(subset[,1] == x), 2]) + cp
  chr.mids <- c(chr.mids, chr.start[y] + cl/2)
  cp = cl + cp + gap
  y <- y + 1
}



setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/00_ITP_bioRxiv_All_Key_Files/11_FiguresDanny")
write.table(lods.f, "BW_females.txt", sep = "\t", quote = FALSE)
write.table(lods.m, "BW_males.txt", sep = "\t", quote = FALSE)
write.table(lods.C, "BW_combined.txt", sep = "\t", quote = FALSE)

pdf("BWmapping.pdf", width = 26, height = 20)
op <- par(mar = c(4,10,2,2))
op <- par(cex = 1.5)

op <- par(mfrow = c(3,1))
plot(c(0, max(chr.start)), y = c(0, 25), t = 'n', ylab = "LOD", xlab = "Chr",xaxt="n", las=2, main = "Body weight Females")
col <- 1
for(y in c("bw42", "bw6", "bw12", "bw18", "bw24")) {
  for(x in c(1:19, "X")) {
    points(subset[which(subset[,1] == x),"cpos"], lods.f[y, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = col)
  }
  col <- col + 1
}
axis(2, at = seq(1, 25, 1), rep("", length(seq(1, 25, 1))))
axis(1, at = chr.mids, paste0(c(1:19,"X")))
legend("topleft", c("bw42", "bw6", "bw12", "bw18", "bw24"), col = 1:5, lwd=2)


plot(c(0, max(chr.start)), y = c(0, 25), t = 'n', ylab = "LOD", xlab = "Chr",xaxt="n", las=2, main = "Body weight Males")
col <- 1
for(y in c("bw42", "bw6", "bw12", "bw18", "bw24")) {
  for(x in c(1:19, "X")) {
    points(subset[which(subset[,1] == x),"cpos"], lods.m[y, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = col)
  }
  col <- col + 1
}
axis(2, at = seq(1, 25, 1), rep("", length(seq(1, 25, 1))))
axis(1, at = chr.mids, paste0(c(1:19,"X")))
legend("topleft", c("bw42", "bw6", "bw12", "bw18", "bw24"), col = 1:5, lwd=2)


plot(c(0, max(chr.start)), y = c(0, 25), t = 'n', ylab = "LOD", xlab = "Chr",xaxt="n", las=2, main = "Body weight Combined")
col <- 1
for(y in c("bw42", "bw6", "bw12", "bw18", "bw24")) {
  for(x in c(1:19, "X")) {
    points(subset[which(subset[,1] == x),"cpos"], lods.C[y, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = col)
  }
  col <- col + 1
}
axis(2, at = seq(1, 25, 1), rep("", length(seq(1, 25, 1))))
axis(1, at = chr.mids, paste0(c(1:19,"X")))
legend("topleft", c("bw42", "bw6", "bw12", "bw18", "bw24"), col = 1:5, lwd=2)
dev.off()


