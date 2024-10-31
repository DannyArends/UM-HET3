setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)

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



off <- 3.65

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
pdf("Figure_1B_RobPref_green.pdf", width = 32, height = 16)
  par(cex=1.6)
  par(cex.axis=1.2)
  op <- par(mar = c(5.1, 4, 4, 5.1))
  plot(x = c(0, cp), y = c(0, 3 + (16 * off)), t = "n", xaxt = "n", yaxt = "n", xlab = "Chromosome", ylab = "")
  yo <- 16*off
  abline(h = yo + 3.65, lty = 2)
  for(x in rev((1:nrow(lods.cM))[c(TRUE, FALSE, FALSE)])){
    i <- 1
    for(chr in c(1:19, "X")){
      onChr <- names(subset[which(subset[,1] == chr), "cpos"])
      xp <- subset[onChr, "cpos"]
      xp <- c(xp[1], xp, xp[length(xp)])
      yp <- c(0, as.numeric(lods.cM[x, onChr]), 0)

      polygon(x = xp, y = yo + yp, 
              col = c(rgb(0, 100, 0, 100, maxColorValue = 255), rgb(0, 100, 0, 100, maxColorValue = 255))[(i %% 2 == 0) + 1], 
              border = "white", lwd=2)
      i <- i + 1
    }
    abline(h = yo, lty = 2)
    yo <- yo - off
  }

  axis(1, at = chr.mids,  c(1:19, "X"))
  axis(2, at = seq(.5 * off, (17 * off), off),  rownames(lods.cM)[(1:nrow(lods.cM))[c(TRUE, FALSE, FALSE)]], las=2)
  axis(3, at = chr.mids,  c(1:19, "X"))
  axis(4, at = seq(.5 * off, (17 * off), off),  rownames(lods.cM)[(1:nrow(lods.cM))[c(TRUE, FALSE, FALSE)]], las=2)
dev.off()


