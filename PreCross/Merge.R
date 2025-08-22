#
# merge.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Merge Sequenom and Monsterplex data
#

setwd("C:/Github/UM-HET3/files")

# Load vcf coded sequenom
sData <- read.table("sequenom/sequenom.data.aligned.txt", sep="\t", row.names=1)
sMap <- read.table("sequenom/sequenom.map.aligned.txt", sep="\t", row.names=1)
sInd <- read.table("sequenom/sequenom.ind.aligned.txt", sep="\t", row.names=1)

# Load vcf coded monsterplex
mData <- read.table("monsterplex/monsterplex.vcfdata.dedup.juli2021.txt", sep="\t", row.names=1)
mMap <- read.table("monsterplex/monsterplex.snp.annot.juli2021.txt", sep="\t", row.names=1)
mInd <- read.table("monsterplex/monsterplex.ind.annot.juli2021.txt", sep="\t", row.names=1)

uniqueInd <- unique(c(rownames(mInd), rownames(sInd)))
uniqueMar <- unique(c(rownames(mMap), rownames(sMap)))

map <- c()
for(x in uniqueMar){
  if(x %in% rownames(sMap)){
    map <- rbind(map, c(x, "Sequenom", sMap[x, 1:2]))
  }else if(x %in% rownames(mMap)){
    map <- rbind(map, c(x, "Monsterplex", mMap[x, 1:2]))
  }else{
    stop()
  }
}
rownames(map) <- map[,1]
map <- map[, -1]
colnames(map) <- c("Origin", "Chr", "Position")

annot <- c()
for(x in uniqueInd){
  if(x %in% rownames(sInd)){
    annot <- rbind(annot, c(x, "Sequenom", sInd[x, c(2)], substr(x,0,2), sInd[x, c(3,4,5,6,7,8,9,10)]))
  }else if(x %in% rownames(mInd)){
    annot <- rbind(annot, c(x, "Monsterplex", NA, mInd[x, 4], mInd[x, c(1,2,3,5,6,7,8,9)]))
  }else{
    stop()
  }
}
rownames(annot) <- annot[,1]
annot <- annot[,-1]
colnames(annot)[c(1,2,3)] <- c("Origin", "Batch", "Site")

merged <- matrix(NA, length(uniqueMar), length(uniqueInd), dimnames=list(uniqueMar, uniqueInd))
for(x in colnames(sData)){
  hasData <- rownames(sData)[which(!is.na(sData[,x]))]
  merged[hasData, x] <- sData[hasData, x]
}

for(x in colnames(mData)){
  hasData <- rownames(mData)[which(!is.na(mData[,x]))]
  #cat(x, "\n")
  merged[hasData, x] <- mData[hasData, x]
}

# Sort by positions
map <- map[sort(as.numeric(map[,"Position"]), index.return=TRUE)$ix, ]
# Sort by the normal chromosome order
smap <- c()
for(chr in c(1:19, "X", "Y", "MT")){
  smap <- rbind(smap, map[which(map[,"Chr"] == chr),])
}
map <- smap

dim(merged)
dim(map)
dim(annot)

merged <- merged[rownames(map), rownames(annot)]

write.table(merged, "merged/all.vcf.sorted.juli2021.txt", quote=FALSE, sep="\t", na = "")
write.table(map, "merged/map.sorted.juli2021.txt", quote=FALSE, sep="\t", na = "")
write.table(annot, "merged/ind.sorted.juli2021.txt", quote=FALSE, sep="\t", na = "")
