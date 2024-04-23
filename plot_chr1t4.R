setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
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

lods.cM <- read.table("progressiveMapping_all.txt", sep = "\t",check.names=FALSE)
lods.mM <- read.table("progressiveMapping_males.txt", sep = "\t",check.names=FALSE)
lods.fM <- read.table("progressiveMapping_females.txt", sep = "\t",check.names=FALSE)
lods.cI <- read.table("progressiveMapping_INT.txt", sep = "\t",check.names=FALSE)

subset <- map[which(map[,1] %in% c(1:19, "X")),]
subset <- cbind(subset, cpos = NA)
gap <- 30000000
chr.start <- c(0)
chr.mids <- c()
cp <- 0
for(chr in c(1:19, "X")){
  cl <- max(as.numeric(subset[which(subset[,1] == chr),2]))
  chr.start <- c(chr.start, cl + cp + gap)
  subset[which(subset[,1] == chr), "cpos"] <- as.numeric(subset[which(subset[,1] == chr), 2]) + cp
  if(chr == "X"){ chr <- 20 }
  cat(chr, " ", cl, " ", cp, " ", gap, " ", chr.start[chr], "\n")
  chr.mids <- c(chr.mids, chr.start[as.numeric(chr)] + cl/2)
  cp = cl + cp + gap
}

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
library(svglite)
pdf("Figure_1_basicQTL.pdf", width = 36, height = 12)
  par(cex=2)
  par(cex.axis=1.5)
  plot(c(0, max(chr.start)), y = c(0, 10), t = 'n', ylab = "LOD", xlab = "Chromosome",xaxt="n", las=2, main = "Longevity (≥ 365 days)")
  for(x in c(1:19, "X")){
    points(subset[which(subset[,1] == x),"cpos"], lods.cM[1, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = "black",lwd=2)
    points(subset[which(subset[,1] == x),"cpos"], lods.mM[1, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = "blue",lwd=2)
    points(subset[which(subset[,1] == x),"cpos"], lods.fM[1, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = "hotpink",lwd=2)
    points(subset[which(subset[,1] == x),"cpos"], lods.cI[1, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = "orange",lwd=2, lty=3)
  }
  abline(h = 3.65, lty=2, col = "green")
  axis(1, at = chr.mids, paste0("", c(1:19, "X")), cex.axis=1.2, las=1)
  legend("topleft", c("All", "Males", "Females", "Interaction"), lwd=2, lty=c(1,1,1,3), col = c("black", "blue", "hotpink", "orange"))
dev.off()



library(svglite)
# Plot the QTL profile
library(RColorBrewer)
colz <- brewer.pal(5, "Purples")

svglite("SeXxMarker_ChrAll.svg", width = 36, height = 12)
par(cex=2)
par(cex.axis=1.5)
layout(matrix(c(1,1,1,1,1,1,2), ncol=7, byrow=TRUE))
op <- par(mar = c(4.5,8,3,0))
image(1:ncol(lods.cI), 1:nrow(lods.cI), t(lods.cI), xlab="4-way map",
      yaxt = 'n',ylab = "",xaxt='n', main="Sex x Marker longevity", breaks = c(0, 2, 3.65, 4.25, 4.95, 5.95, 100), col=c("white", colz))

abline(v=which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
yaxisD <- paste0(">", seq(365, 1100, 15), " days")
axis(2, at = (1:length(yaxisD))[seq(1,length(yaxisD), 4)], yaxisD[seq(1,length(yaxisD), 4)], las = 2)
chrP <- c(0, which(diff(as.numeric(as.factor(map[,"Chr"]))) != 0))
axis(1, at = chrP + (diff(c(chrP, nrow(map)))/2), c(1:19, "X"))
box()
op <- par(mar = c(4,0,3,0))
plot(1:10, t = 'n', bty="n", xaxt ='n', yaxt ='n',xlab="")
legend("top", legend= c("<2.00", "2.00 - 3.65", "3.65 - 4.25", "4.25 - 4.95", "4.95 - 5.95", ">5.95"), fill = c("white", colz), bg="white")

dev.off()


