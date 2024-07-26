timepoints c(365, 560, 875, 1025)
#Vita2c (females)
#Vita9b (males)

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
tp <- "> 365"
#timepoints c(365, 560, 875, 1025)

pdf("Figure_1B_TPS_365.pdf", width = 36, height = 12)
  #par(cex=2)
  #par(cex.axis=1.5)
  plot(c(1497316549 + 10000000, 1651373135 - 30000000), y = c(0, 6), t = 'n', ylab = "LOD", xlab = "Chromosome",xaxt="n", las=2, main = "Longevity (365 days)")
  for(x in c("9")){
    points(subset[which(subset[,1] == x),"cpos"], lods.cM[tp, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = "black",lwd=2)
    points(subset[which(subset[,1] == x),"cpos"], lods.mM[tp, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = "blue",lwd=2)
    points(subset[which(subset[,1] == x),"cpos"], lods.fM[tp, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = "hotpink",lwd=2)
    points(subset[which(subset[,1] == x),"cpos"], lods.cI[tp, rownames(subset)[which(subset[,1] == x)]], t = 'l', col = "orange",lwd=2, lty=3)
  }
  abline(h = 3.65, lty=2, col = "green")
  axis(1, at = chr.mids, paste0("", c(1:19, "X")), cex.axis=1.2, las=1)
  legend("topleft", c("All", "Males", "Females", "Interaction"), lwd=2, lty=c(1,1,1,3), col = c("black", "blue", "hotpink", "orange"))
dev.off()

