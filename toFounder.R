setwd("C:/Github/UM-HET3/files/merged")

# Load vcf coded merged data
gtC <- read.table("all.vcf.sorted.txt", sep="\t", row.names=1,colClasses="character")
map <- read.table("map.sorted.txt", sep="\t", row.names=1)
ind <- read.table("ind.sorted.txt", sep="\t", row.names=1)

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

write.table(gtC, "selected.vcf.sorted.txt", quote=FALSE, sep="\t", na = "")

gts.cor <- cor(t(gts), use="pair")
chr1.names <- rownames(map)[which(map[,"Chr"] == 1 & map[,"Origin"] == "Monsterplex")]
image(gts.cor[chr1.names,chr1.names])
chr1.names <- rownames(map)[which(map[,"Chr"] == 1 & map[,"Origin"] == "Sequenom")]
image(gts.cor[chr1.names,chr1.names])

# Read the founder snps
fvcf <- read.table("fvcfAll.txt", sep="\t", row.names=1,colClasses="character")

# Bind the position names (since rsIDs are not in the founders), and test if in founders
map <- cbind(map, posname = paste0(map[, "Chr"], "_", map[, "Position"]))
map <- cbind(map, inFounder = as.numeric(map[, "posname"] %in% rownames(fvcf)))
map <- cbind(map, type = NA)

# 4way autosomes, X chromosome and matrix
m4wayA <- mok[map[, "inFounder"] == 1 & map[, "Chr"] != "X"]
m4way <- m4wayA

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

table(unlist(gts4way))

tbls.mar <- apply(gts4way, 1, table)
tbls.ind <- apply(gts4way, 2, table)

wrong <- names(which(unlist(lapply(tbls.mar, function(x){ as.numeric(x["XX"]) / sum(x) })) > 0.1))
wrong <- c(wrong, names(which(mar.missing(gts4way) < 0.1)))

mOK <- rownames(gts4way)[which(!rownames(gts4way) %in% wrong)]

gts4way <- gts4way[mOK,]
map <- map[mOK, ]

write.table(gts4way, "gts4way.txt", quote=FALSE, sep="\t")
write.table(map, "map.gts4way.txt", quote=FALSE, sep="\t")
write.table(ind, "ind.gts4way.txt", quote=FALSE, sep="\t")

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

write.table(gts4way, "gts4way.rqtl.txt", quote=FALSE, sep="\t", na = "")
