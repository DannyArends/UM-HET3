#
# PostPhasingRemoval.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Post phasing removal of markers not in linkage with the neighbouring markers
#

setwd("C:/Github/UM-HET3/files")
gts4way.fill <- read.table("merged/gts4way.rqtl.filled.Juli2021.txt", sep="\t")
gts4way.rqtl <- read.table("merged/gts4way.rqtl.Juli2021.txt", sep="\t")

# Use correlation to catch weird markers
mcor <- cor(t(gts4way.fill[,-c(1:2)]), use="pair")
op <- par(mfrow=c(1, 2))
image(1:nrow(mcor), 1:ncol(mcor), mcor, xlab="Marker", ylab="Marker", main="Correlation (After Imputation) All markers")
badM <- names(which(apply(apply(mcor,1,is.na),1,sum) > 50)) # Some NAs are fine, but more than 50 means you're a bad marker

gts4way.rqtl <- gts4way.rqtl[-which(rownames(gts4way.rqtl) %in% badM),]
gts4way.fill <- gts4way.fill[-which(rownames(gts4way.fill) %in% badM),]

mcor <- cor(t(gts4way.fill[,-c(1:2)]), use="pair")
image(1:nrow(mcor), 1:ncol(mcor), mcor, xlab="Marker", ylab="Marker", main="Correlation (After Imputation) All markers")

# Zoom into a chromosome
chr <- 1
image(which(gts4way.fill[,1] == chr), which(gts4way.fill[,1] == chr), mcor[gts4way.fill[,1] == chr, gts4way.fill[,1] == chr])
box()

# Chromsome 6 weird 434:437
# Chromsome 7 weird 494:497
# Chromsome 12 weird 849:850
# Chromsome 17 weird 1061:1065

ii <- c(429)
gts4way.rqtl <- gts4way.rqtl[-ii,]
gts4way.fill <- gts4way.fill[-ii,]

# too much missing data
tooMM <- names(which(apply(apply(gts4way.fill[-c(1,2),],1,is.na),2,sum)/ncol(gts4way.fill) > 0.7))

gts4way.rqtl <- gts4way.rqtl[-which(rownames(gts4way.rqtl) %in% tooMM),]
gts4way.fill <- gts4way.fill[-which(rownames(gts4way.fill) %in% tooMM),]

mcor <- cor(t(gts4way.fill[,-c(1:2)]), use="pair")
image(1:nrow(mcor), 1:ncol(mcor), mcor, xlab="Marker", ylab="Marker", main="Correlation (After Imputation) All markers, Cleaned")

write.table(gts4way.fill, "merged/gts4way.filled.cleaned.GN2.txt", quote=FALSE, sep="\t", na = "")
write.table(gts4way.rqtl, "merged/gts4way.rqtl.filled.cleaned.txt", quote=FALSE, sep="\t", na = "")

map <- read.table("merged/map.gts4way.juli2021.txt", sep="\t")
map <- map[rownames(gts4way.rqtl),]
write.table(map, "merged/map.gts4way.filled.cleaned.GN2.txt", quote=FALSE, sep="\t", na = "")

umhet3.geno <- cbind(Chr = map[,"Chr"], Locus = rownames(map), cM = NA, Mb = map[, "Position"] / 1000000, gts4way.rqtl)
write.table(umhet3.geno, "UMHET3.Juli21.geno", quote=FALSE, sep="\t", na = "",row.names=FALSE)

# No duplicated markers (merge sequenom with monsterplex individuals
umhet3.geno.unique <- NULL
x <- 1
while(x < nrow(map)){
  cat("Marker: ", x, "\n")
  ii <- which(map[,"posname"] == map[x,"posname"])
  if(length(ii) == 1){
    umhet3.geno.unique <- rbind(umhet3.geno.unique, umhet3.geno[x, ])
  }else{
    cat("Iz: ", ii, "\n")
    gts <- apply(umhet3.geno[ii,-c(1:4)],2,function(x){
      if(all(is.na(x))){return(NA)}
      if(length(which(is.na(x))) == 1) return(na.omit(x))
      if(length(which(is.na(x))) == 0 && x[1] == x[2]) return(x[1])
      if(length(which(is.na(x))) == 0 && x[1] != x[2]){ cat("!"); return(NA) }
    })
    umhet3.geno.unique <- rbind(umhet3.geno.unique, c(umhet3.geno[ii[1],c(1:4)], gts))
    x <- max(ii)
  }
  x <- x + 1
}
rownames(umhet3.geno.unique) <- umhet3.geno.unique[, "Locus"]
write.table(umhet3.geno.unique, "UMHET3.Juli21.noDup.geno", quote=FALSE, sep="\t", na = "",row.names=FALSE)

# Create cross object for determining genotypes
setwd("C:/Github/UM-HET3/files")
ind <- read.table("merged/ind.gts4way.Juli2021.txt", sep="\t")

write.table(cbind(NA,NA, t(cbind(Individual = rownames(ind), ind))), file = "UMHET3.Juli21.noDup.geno", col.names = FALSE, sep = ",", quote=FALSE, na = "")
write.table(rbind(GenoID = c(NA, NA, colnames(umhet3.geno.unique)[-c(1:4)]), 
                  cbind(Chr = umhet3.geno.unique$Chr, Mb = as.numeric(umhet3.geno.unique$Mb) , umhet3.geno.unique[, -c(1:4)])
                 ), file = "UMHET3.Juli21.noDup.geno", col.names = FALSE, sep = ",", append=TRUE, quote=FALSE, na="")

# Read cross object for determining genotypes
library(qtl)
mcross <- read.cross(format="csvr", file="UMHET3.Juli21.noDup.geno", genotypes=NULL, na.strings=c("-", "NA"))
nxo <- countXO(mcross)
plot(nxo)


setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/merged/")
ind <- read.table("ind.gts4way.Juli2021.txt", sep="\t")
umhet3.geno.unique <- read.table("UMHET3.Juli21.noDup.geno", sep=",")

write.table(cbind(NA,NA, t(cbind(Individual = rownames(ind), ind))), file = "UMHET3.Juli21.noDup.1.geno", col.names = FALSE, sep = ",", quote=FALSE, na = "")
write.table(rbind(GenoID = c(NA, NA, colnames(umhet3.geno.unique)[-c(1:4)]), 
                  cbind(Chr = umhet3.geno.unique$Chr, Mb = as.numeric(umhet3.geno.unique$Mb) , umhet3.geno.unique[, -c(1:4)])
                 ), file = "UMHET3.Juli21.noDup.1.geno", col.names = FALSE, sep = ",", append=TRUE, quote=FALSE, na="")



