#
# monsterplex.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Merge Sequenom genotype data with ITP phenotype and covariate data, order the map, and transform genotypes to vcf coding (similar to monsterplex data)
#

setwd("C:/Github/UM-HET3/files/sequenom")
sdata <- read.csv("sequenom.data.csv", row.names=1, na.strings=c("U", "x", "", "NA"),colClasses="character")
aind <- read.csv("sequenom.ind.annot.txt",sep="\t", na.strings=c("U", "x", "", "NA"))
amar <- read.csv("sequenom.snp.annot.txt",sep="\t")

# Some individuals do not have any data, remove them
aind <- aind[-which(is.na(aind[,4])),]
rownames(aind) <- aind[,4]
aind <- aind[, -c(1,2,5)]

# Remove the individuals for whom we have no annotation
sdata <- sdata[-which(!rownames(sdata) %in% rownames(aind)),]

# Sort the data by the annotation, and use the UM-HET3 ID
sdata <- sdata[rownames(aind),]
rownames(sdata) <- aind[,1]
rownames(aind) <- aind[,1]
aind <- aind[,-1]

map <- c()
for(x in colnames(sdata)){
  pos <- amar[which(amar[,1] == x)[1],2]
  pos <- strsplit(pos, "-")[[1]][1]
  map <- rbind(map, c(x, strsplit(pos, ":")[[1]][1], strsplit(pos, ":")[[1]][2], amar[which(amar[,1] == x)[1],3]))
}

# Sort by positions
map <- map[sort(as.numeric(map[,3]), index.return=TRUE)$ix, ]
# Sort by the normal chromosome order
smap <- c()
for(chr in c(1:19, "X", "Y", "MT")){
  smap <- rbind(smap, map[which(map[,2] == chr),])
}
rownames(smap) <- smap[,1]
smap <- smap[,-1]
colnames(smap) <- c("Chr", "Position", "Allele")

# Order genotypes as normal
sdata <- sdata[, rownames(smap)]
# Code the genotypes as normal
sdata[sdata == "A"] <- "AA"
sdata[sdata == "C"] <- "CC"
sdata[sdata == "G"] <- "GG"
sdata[sdata == "T"] <- "TT"
sdata <- t(sdata)

vcfCoded <- matrix(NA, nrow(sdata), ncol(sdata), dimnames=list(rownames(sdata), colnames(sdata)))
nonS <- c()
for(x in rownames(smap)){
  alt <- smap[x, "Allele"]
  mtable <- table(unlist(strsplit(sdata[x,], "")))
  if(length(mtable) == 1){
    nonS <- c(nonS, x)
  }else{
    ref <- names(mtable)[which(!names(mtable) %in% alt)]
    if(length(ref) == 1){
      cat(ref, alt, "\n")
      isRef <- which(sdata[x,] == paste0(ref, ref))
      isAlt <- which(sdata[x,] == paste0(alt, alt))
      isHet <- which(sdata[x,] == paste0(ref, alt) | sdata[x,] == paste0(alt, ref))
      vcfCoded[x, isRef] <- "0/0"
      vcfCoded[x, isHet] <- "0/1"
      vcfCoded[x, isAlt] <- "1/1"
    }else{
      nonS <- c(nonS, x)
    }
  }
}
vcfCoded <- vcfCoded[-which(rownames(vcfCoded) %in% nonS),]
smap <- smap[rownames(vcfCoded),]

write.table(vcfCoded, "sequenom.data.aligned.txt", quote=FALSE, sep="\t")
write.table(smap, "sequenom.map.aligned.txt", quote=FALSE, sep="\t")
write.table(aind, "sequenom.ind.aligned.txt", quote=FALSE, sep="\t")
