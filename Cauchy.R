cauchycombination <- function(p) {
   qc <- qcauchy(p, lower.tail=FALSE)
   pc <- pcauchy(mean(qc), lower.tail=FALSE)
   pc
}

setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

lods <- read.table("progressiveMapping_all.txt", sep = "\t", check.names = FALSE)

cc <- c()
for(x in 1:ncol(lods)){
  cc <- c(cc, cauchycombination(10 ^ -as.numeric(lods[,x])))
}

names(cc) <- colnames(lods)

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

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

map <- cbind(Chr = chrs, Pos = positions)
rownames(map) <- colnames(pull.geno(mcross))

plot(c(0, max(chr.start)), y = c(0, 10), t = 'n', ylab = "LOD", xlab = "Chromosome",xaxt="n", las=2, main = "Longevity (≥ 365 days)")
for(x in c(1:19, "X")){
  points(subset[which(subset[,1] == x),"cpos"], -log10(cc[rownames(subset)[which(subset[,1] == x)]]), t = 'l', col = "black",lwd=2)
}
abline(h = 3.65, lty=2, col = "green")
axis(1, at = chr.mids, paste0("", c(1:19, "X")), cex.axis=1.2, las=1)
legend("topleft", c("All", "Males", "Females", "Interaction"), lwd=2, lty=c(1,1,1,3), col = c("black", "blue", "hotpink", "orange"))

all <- c("1_3010272", "1_24042124", "1_120474787", "2_89844287", "2_112712327", "2_139956785","3_83838529", "3_92135706", "4_52524395",
         "5_67573068", "6_107382038", "6_140010438", "9_29939029", "9_104091597", "9_124056586", "10_72780332", "11_5628810", "11_82178599",
         "12_112855820", "13_89689878", "14_101437457", "15_74248242", "17_32883804", "18_60822951", "X_36008085", "X_156343080")


cat(p.adjust(cc, "BH")[all], sep = "\n")

