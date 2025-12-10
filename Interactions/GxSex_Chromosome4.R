#
# GxS_chr4.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Code to produce a visualization of G x Sex effects on chromosome 4 (figure not used in main or supplements)
#

library(qtl)

source("ActuarialMapping/adjustXprobs.R")

# Read cross object
mcross <- read.cross(format="csvr", file="DataSet/um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

# Writing out 
iix <- which(pull.pheno(mcross)[,"longevity"] > 365)
write.table(pull.pheno(mcross)[iix,"GenoID"], "Cases_UM_HET3.txt",sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)

# Create the map object
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

c4 <- which(map[,1] == "4")

setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
lods.i <- read.table("DataSet/output/progressiveMapping_INT.txt", sep = "\t", check.names = FALSE)
lods.cM <- read.table("DataSet/output/progressiveMapping_all.txt", sep = "\t",check.names=FALSE)

pdf("DataSet/output/GxSonChr4.pdf", width = 24, height = 12)
 op <- par(mar = c(10,5,2,2))
 par(cex=1.5)
 image(x = 1:length(c4), y = 1:nrow(lods.i), t(as.matrix(lods.i[, c4])), xaxt="n", yaxt="n", xlab = "", ylab = "time", main = "GxS interaction")
 axis(1, at = 1:length(c4), rownames(map[c4,]), las=2)
 axis(2, at = 1:nrow(lods.i), seq(365,1100,15), las=2)
dev.off()

pdf("DataSet/output/Chr4_Vita4aorTwo.pdf", width = 24, height = 12)
 op <- par(mar = c(10,5,2,2))
 par(cex=1.5)
 image(x = 1:length(c4), y = 1:nrow(lods.cM), t(as.matrix(lods.cM[, c4])), xaxt="n", yaxt="n", xlab = "", ylab = "time", main = "Main effect / No interaction")
 axis(1, at = 1:length(c4), rownames(map[c4,]), las=2)
 axis(2, at = 1:nrow(lods.cM), seq(365,1100,15), las=2)
dev.off()

