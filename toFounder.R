setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/merged")

# Load vcf coded merged data
gtC <- read.table("all.vcf.sorted.juli2021.txt", sep="\t", row.names=1,colClasses="character")
map <- read.table("map.sorted.juli2021.txt", sep="\t", row.names=1)
ind <- read.table("ind.sorted.juli2021.txt", sep="\t", row.names=1)

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/")
cat(unique(rownames(ind)), file = "Genotypes_Total.txt", sep = "\n")

gts <- apply(gtC,2,function(x){ return(as.numeric(factor(x, levels=c("0/0", "0/1", "1/1")))) } )
rownames(gts) <- rownames(gtC)

mar.missing <- function(gts){ apply(gts, 1, function(x){ length(which(!is.na(x))) / length(x) }) }
ind.missing <- function(gts){ apply(gts, 2, function(x){ length(which(!is.na(x))) / length(x) }) }

mok <- names(which(mar.missing(gts) > 0.1))
iok <- names(which(ind.missing(gts[mok,]) > 0.1))
mok <- names(which(mar.missing(gts[mok,iok]) > 0.1))

cat(length(mok), "Markers", length(iok), "individuals\n")
gts <- gts[mok, iok]
gtC <- gtC[mok, iok]
map <- map[mok,]
ind <- ind[iok,]

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/")
cat(unique(rownames(ind)), file = "Max_90_Perc.txt", sep = "\n")


write.table(gtC, "selected.vcf.sorted.juli2021.txt", quote=FALSE, sep="\t", na = "")

gts.cor <- cor(t(gts), use="pair")
chr1.names <- rownames(map)[which(map[,"Chr"] == 1 & map[,"Origin"] == "Monsterplex")]
image(gts.cor[chr1.names,chr1.names])
chr1.names <- rownames(map)[which(map[,"Chr"] == 1 & map[,"Origin"] == "Sequenom")]
image(gts.cor[chr1.names,chr1.names])

# Read the founder snps
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files/merged")
fvcf <- read.table("fvcfAll.txt", sep="\t",colClasses="character", header=TRUE)
rownames(fvcf) <- paste0(fvcf[,1], "_", fvcf[,2])
#setwd("C:/Github/UM-HET3/files/merged")

# Bind the position names (since rsIDs are not in the founders), and test if in founders
map <- cbind(map, posname = paste0(map[, "Chr"], "_", map[, "Position"]))
map <- cbind(map, inFounder = as.numeric(map[, "posname"] %in% rownames(fvcf)))
map <- cbind(map, type = NA)

# 4way autosomes, X chromosome and matrix
m4wayA <- mok[map[, "inFounder"] == 1 & map[, "Chr"] != "X"]
m4wayX <- mok[map[, "inFounder"] == 1 & map[, "Chr"] == "X"]
m4way <- c(m4wayA,m4wayX)

gts4way <- matrix(NA, length(m4way), ncol(gtC), dimnames= list(m4way, colnames(gtC)))
cat(length(m4way), "Markers", ncol(gtC), "individuals\n")

founders <- colnames(fvcf)[6:9]

for (m in m4wayA) {
  pm <- map[m, "posname"] # position name 
  if(sum(fvcf[pm, founders] == "1/1") == 1){ # One founder
    if (fvcf[pm, "BALB_cJ"] == "1/1") {
      gts4way[m, which(gtC[m,] == "0/0")] <-  "B?"
      gts4way[m, which(gtC[m,] == "0/1")] <-  "A?"
      gts4way[m, which(gtC[m,] == "1/1")] <-  "XX"
      map[m, "type"] = "Mat1"
    }
    if (fvcf[pm, "C3H_HeJ"] == "1/1") {
      gts4way[m, which(gtC[m,] == "0/0")] <-  "?D"
      gts4way[m, which(gtC[m,] == "0/1")] <-  "?C"
      gts4way[m, which(gtC[m,] == "1/1")] <-  "XX"
      map[m, "type"] = "Pat1"
    }
    if(fvcf[pm, "DBA_2J"] == "1/1") {
      gts4way[m, which(gtC[m,] == "0/0")] <-  "?C"
      gts4way[m, which(gtC[m,] == "0/1")] <-  "?D"
      gts4way[m, which(gtC[m,] == "1/1")] <-  "XX"
      map[m, "type"] = "Pat2"
    }
  }
  if(sum(fvcf[pm, founders] == "1/1") == 2){ # 2 Founders
    if(fvcf[pm, "BALB_cJ"] == "1/1") {
      if(fvcf[pm, "C3H_HeJ"] == "1/1"){
        gts4way[m, which(gtC[m,] == "0/0")] <-  "BD"
        gts4way[m, which(gtC[m,] == "0/1")] <-  "??"
        gts4way[m, which(gtC[m,] == "1/1")] <-  "AC"
        map[m, "type"] = "Pat1"
      }
      if(fvcf[pm, "DBA_2J"] == "1/1"){
        gts4way[m, which(gtC[m,] == "0/0")] <-  "BC"
        gts4way[m, which(gtC[m,] == "0/1")] <-  "??"
        gts4way[m, which(gtC[m,] == "1/1")] <-  "AD"
        map[m, "type"] = "Pat2"
      }
    }else{
      cat("Marker", m, "cannot be decomposed into founder genotypes\n")
      map[m, "type"] = "U"
    }
  }
  if(sum(fvcf[pm, founders] == "1/1") == 3){ # All founders (except for the B6 have a SNP)
    gts4way[m, which(gtC[m,] == "1/1")] <-  "A?"
    gts4way[m, which(gtC[m,] == "0/1")] <-  "B?"
    gts4way[m, which(gtC[m,] == "0/0")] <-  "XX"
    map[m, "type"] = "Mat2"
  }
  cat("Done", m, "=", pm, "\n")
}

males <- ind[colnames(gtC),"Sex"] == "M"
females <- ind[colnames(gtC),"Sex"] == "F"

for (m in m4wayX) { # Only two posibilities are observed in the founder conversion map
  pm <- map[m, "posname"] # position name
  if(sum(fvcf[pm, founders] == "1/1") == 1){ # One founder
    if (fvcf[pm, "BALB_cJ"] == "1/1") {
      gts4way[m, which(gtC[m,] == "0/0" & males)] <-  "BD"
      gts4way[m, which(gtC[m,] == "0/1" & males)] <-  "AD"
      gts4way[m, which(gtC[m,] == "1/1" & males)] <-  "AD"
      gts4way[m, which(gtC[m,] == "0/0" & females)] <-  "BC"
      gts4way[m, which(gtC[m,] == "0/1" & females)] <-  "AC"
      gts4way[m, which(gtC[m,] == "1/1" & females)] <-  "XX"
    }
  }else if(sum(fvcf[pm, founders] == "1/1") == 3){ # 3 founders
      gts4way[m, which(gtC[m,] == "0/0" & males)] <-  "BD"
      gts4way[m, which(gtC[m,] == "0/1" & males)] <-  "XX"
      gts4way[m, which(gtC[m,] == "1/1" & males)] <-  "AD"
      
      gts4way[m, which(gtC[m,] == "0/0" & females)] <-  "XX"
      gts4way[m, which(gtC[m,] == "0/1" & females)] <-  "BC"
      gts4way[m, which(gtC[m,] == "1/1" & females)] <-  "AC"
  }
  cat("Done", m, "=", pm, "\n")
}

table(unlist(gts4way))

tbls.mar <- apply(gts4way, 1, table)
tbls.ind <- apply(gts4way, 2, table)

wrong <- names(which(unlist(lapply(tbls.mar, function(x){ as.numeric(x["XX"]) / sum(x) })) > 0.1))
wrong <- c(wrong, names(which(mar.missing(gts4way) < 0.1)))

mOK <- rownames(gts4way)[which(!rownames(gts4way) %in% wrong)]

gts4way <- gts4way[mOK,]
map <- map[mOK, ]

write.table(gts4way, "gts4way.Juli2021.txt", quote=FALSE, sep="\t",na="")
write.table(map, "map.gts4way.Juli2021.txt", quote=FALSE, sep="\t",na="")
write.table(ind, "ind.gts4way.Juli2021.txt", quote=FALSE, sep="\t",na="")

gts4way[gts4way=="AC"] <- 1
gts4way[gts4way=="BC"] <- 2
gts4way[gts4way=="AD"] <- 3
gts4way[gts4way=="BD"] <- 4
gts4way[gts4way=="A?"] <- 5
gts4way[gts4way=="A-"] <- 5
gts4way[gts4way=="B?"] <- 6
gts4way[gts4way=="B-"] <- 6
gts4way[gts4way=="?C"] <- 7
gts4way[gts4way=="?D"] <- 8
gts4way[gts4way=="XX"] <- NA
gts4way[gts4way=="??"] <- NA

gts4num <- apply(gts4way, 2, as.numeric)
rownames(gts4num) <- rownames(gts4way)
table(unlist(gts4way))

#write.table(gts4way, "gts4way.rqtl.Juli2021.txt", quote=FALSE, sep="\t", na = "")

gts4wayGN2 <- cbind(map[rownames(gts4way),c(2,3)], gts4way)
write.table(gts4wayGN2, "gts4way.rqtl.GN2.Juli2021.txt", quote=FALSE, sep="\t", na = "")


gts4wayGN2 <- cbind(posname = paste0(gts4wayGN2[,"Chr"], "_", gts4wayGN2[, "Position"]), gts4wayGN2)

# No duplicated markers (merge sequenom rs markers with monsterplex sequencing markers)
umhet3.geno.unique <- NULL
x <- 1
while(x < nrow(gts4wayGN2)){
  cat("Marker: ", x, "\n")
  ii <- which(gts4wayGN2[,"posname"] == gts4wayGN2[x,"posname"])
  if(length(ii) == 1){
    umhet3.geno.unique <- rbind(umhet3.geno.unique, gts4wayGN2[x, ])
  }else{
    cat("Iz: ", ii, "\n")
    gts <- apply(gts4wayGN2[ii,-c(1:3)],2,function(x){
      if(all(is.na(x))){return(NA)}
      if(length(which(is.na(x))) == 1) return(na.omit(x))
      if(length(which(is.na(x))) == 0 && x[1] == x[2]) return(x[1])
      if(length(which(is.na(x))) == 0 && x[1] != x[2]){ cat("!"); return(NA) }
    })
    umhet3.geno.unique <- rbind(umhet3.geno.unique, c(gts4wayGN2[ii[1],c(1:3)], gts))
    x <- max(ii)
  }
  x <- x + 1
}
rownames(umhet3.geno.unique) <- umhet3.geno.unique[,1]
umhet3.geno.unique <- umhet3.geno.unique[,-1]
write.table(umhet3.geno.unique, "gts4way.rqtl.GN2.Juli2021.NoDUP.txt", quote=FALSE, sep="\t", na = "")

# Create cross object for determining genotypes
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/merged")
write.table(cbind(NA,NA, t(cbind(Individual = rownames(ind), ind))), file = "um-het3-rqtl.noDup.csvr", col.names = FALSE, sep = ",", quote=FALSE, na = "")
write.table(rbind(GenoID = c(NA, NA, colnames(umhet3.geno.unique)[-c(1,2)]), 
                  cbind(Chr = umhet3.geno.unique[,"Chr"], Mb = as.numeric(umhet3.geno.unique[,"Position"]) / 1000000, umhet3.geno.unique[, -c(1,2)])
                 ), file = "um-het3-rqtl.noDup.csvr", col.names = FALSE, sep = ",", append=TRUE, quote=FALSE, na="")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.noDup.csvr", genotypes=NULL, na.strings=c("-", "NA"))

# Compute crossovers
nxo <- countXO(mcross)
plot(nxo)

# Subset the cross using only the individuals with less than 40 XOs
goodInd <- which(nxo <= 80)
gcross <- subset(mcross, ind = goodInd)

ind <- pull.pheno(gcross)[,1]
geno <- pull.geno(gcross)
rownames(geno) <- ind
gfull <- cbind(umhet3.geno.unique[colnames(geno), 1:2], t(geno))

badMarker <- "4_58784823" # Remove the final misbehaving marker on chromosome 4
gfull <- gfull[-which(rownames(gfull) == badMarker),]
dim(gfull)

write.table(gfull, file="ITP_6480x893_MatPat_Sep21.txt", sep = "\t", quote=FALSE,na="")

# Paternal and Maternal maps
mMap <- unique(map[grep("Mat", map[,"type"]), "posname"])
fMap <- unique(map[grep("Pat", map[,"type"]), "posname"])

gMat <- gfull[which(rownames(gfull) %in% mMap),]
gPat <- gfull[which(rownames(gfull) %in% fMap),]

write.table(gMat, file="ITP_6480x486_Maternal_Sep21.txt", sep = "\t", quote=FALSE,na="")
write.table(gPat, file="ITP_6480x396_Paternal_Sep21.txt", sep = "\t", quote=FALSE,na="")


### Some checks on the new maps, the code below can be ignored

# Quick n dirty QTL scan, using site and year on both data sets
addcovar <- cbind(as.numeric(pull.pheno(mcross)[, "Sex"]), as.numeric(pull.pheno(mcross)[, "Site"]), as.numeric(pull.pheno(mcross)[, "Cohort.Year"]))
allI <- scanone(mcross, pheno.col = "Longevity_HET3_ITP", addcovar=addcovar, method="hk")

addcovar <- cbind(as.numeric(pull.pheno(gcross)[, "Sex"]), as.numeric(pull.pheno(gcross)[, "Site"]), as.numeric(pull.pheno(gcross)[, "Cohort.Year"]))
lessI <- scanone(gcross, pheno.col = "Longevity_HET3_ITP", addcovar=addcovar, method="hk")

# No major differences are seen
plot(allI, lessI)


# Quick and Dirty Permutations
addcovar <- cbind(as.numeric(pull.pheno(gcross)[, "Site"]), as.numeric(pull.pheno(gcross)[, "Cohort.Year"]))
perm.res <- scanone(gcross, pheno.col = "Longevity_HET3_ITP", addcovar=addcovar, method="hk", n.perm=1000, n.cluster=3)

# Find old mice with too many xo
tooManyXO <- pull.pheno(mcross)[which(nxo > 40),]
oldHighXO <- which(tooManyXO[, "Longevity_HET3_ITP"] > 900)
length(oldHighXO)
write.table(pull.pheno(mcross)[which(nxo > 40), c(1, 4,5,6)], file="toReGeno.txt",sep="\t", row.names=FALSE, quote=FALSE)


# Fill genotypes
fcross <- fill.geno(mcross, method = "maxmarginal", error.prob = 0.01, min.prob = 0.85)
gts4way.full <- t(pull.geno(fcross))
colnames(gts4way.full) <- pull.pheno(fcross)[,"Individual"]
# Write out the filled genotypes
gts4way.full <- cbind(map[rownames(gts4way.full),c(2,3)], gts4way.full)
write.table(gts4way.full, "merged/gts4way.rqtl.filled.Juli2021.txt", quote=FALSE, sep="\t", na="")
