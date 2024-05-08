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



cc["1_3010272"]
cc["1_24042124"]
cc["1_132295971"]
cc["2_112712327"]
cc["2_139956785"]
cc["2_157112564"]
cc["3_83838529"]
cc["3_92135706"]
cc["4_55012301"]
cc["4_74811205"]
cc["4_145301445"]
cc["5_67573068"]
cc["6_108075853"]
cc["8_36994142"]
cc["9_29939029"]
cc["9_54904313"]
cc["9_104091597"]
cc["9_124029281"]
cc["10_72780332"]
cc["11_5628810"]
cc["11_82178599"]
cc["12_112855820"]
cc["13_89689878"]
cc["14_101437457"]
cc["15_74248242"]
cc["17_32883804"]
cc["17_18001459"]
cc["17_68770703"]
cc["18_52488251"]
cc["X_36008085"]
cc["X_156343080"]


