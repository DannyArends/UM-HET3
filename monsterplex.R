setwd("C:/Github/UM-HET3/files/monsterplex")

# Load the sample annot
annot <- read.csv("Sample_Metadata.csv", sep=",", na.strings=c("blank", "U", "x", "", "NA"))
annot <- annot[which(annot[,3] == "Monsterplex"),]

# Load the data (children and founders)
cvcf <- read.csv("cvcfAll.txt", sep="\t", na.strings=c("U", "x", "", "NA"))
fvcf <- read.csv("fvcfAll.txt",sep="\t", na.strings=c("U", "x", "", "NA"))

# Match between them
rownames(cvcf) <- paste0(cvcf[,1], "_", cvcf[,2])
cvcf <- cvcf[which(rownames(cvcf) %in% rownames(fvcf)),]

# Remove the annotation from the data file
map <- cvcf[,1:4]
colnames(map)[1:2] <- c("Chr", "Position")
cvcf <- cvcf[,-c(1:5)]

# Remove the name part we don't need
colnames(cvcf) <- gsub(".fastq.gz_default.sorted.bam", "", colnames(cvcf))
colnames(cvcf) <- gsub(".fq.gz_default.sorted.bam", "", colnames(cvcf))
colnames(cvcf) <- gsub(".fastq_default.sorted.bam", "", colnames(cvcf))
colnames(cvcf) <- gsub("..Mouse_", "", colnames(cvcf))
colnames(cvcf) <- gsub("..", "", colnames(cvcf),fixed=TRUE)

# Split out the UM and JL samples (their UM-Het3 ID is in the name)
UMs <- grep("UM", colnames(cvcf))
JLs <- grep("JL", colnames(cvcf))
noA <- 1:ncol(cvcf)
noA <- noA[-which(noA %in% c(UMs, JLs))]

# Find out what the samples are that remain, and build mAnnot
mAnnot <- cbind(colnames(cvcf)[noA], UM_HET3ID = NA)
cnt <- 1
for(id in mAnnot[,1]){
  idx <- strsplit(id, "_")[[1]]
  #cat(paste0(idx, collapse=" "), "\n")
  if(id == "M1toM5_Um14065"){
    mAnnot[cnt,2] <- "UM14065"
  }else if(length(idx) == 5){
    mID <- idx[1]
    sID <- idx[3]
    inAnnot <- which(grepl(sID, annot[, "Source"]) & grepl(mID, annot[, "Source"]))
    if(length(inAnnot) == 1){
      mAnnot[cnt,2] <- annot[inAnnot,"ITP_ID"]
    }else{
      pGenoID <- paste0(idx[3],"_",idx[4],"_",idx[5])
      inAnnot <- which(annot[,"Geno_ID"] == pGenoID)
      if(length(inAnnot) == 1){
        mAnnot[cnt,2] <- annot[inAnnot,"ITP_ID"]
      }else{
        stop()
      }
    }
  }else if(length(idx) == 4){
    pGenoID <- paste0(idx[3], "_", idx[4])
    inAnnot <- which(annot[,"Geno_ID"] == pGenoID)
    if(length(inAnnot) == 1){
      mAnnot[cnt,2] <- annot[inAnnot,"ITP_ID"]
    }else{
      stop()
    }
  }
  cnt <- cnt + 1
}

# Add the UM samples
UMids <- colnames(cvcf[, UMs])
UMidx <- gsub("M1toM5_", "", UMids)
UMidx <- gsub("Extra_40_markers_M1toM12_OUTPUT_", "", UMidx)
mAnnot <- rbind(mAnnot, cbind(UMids, UMidx))

# Add the JL samples
JLids <- colnames(cvcf[, JLs])
JLidx <- gsub("M21to24data_", "", JLids)
JLidx <- gsub("M12to20_Aug3_", "", JLidx)
JLidx <- gsub("M11_M15redone_", "", JLidx)
JLidx <- gsub("Extra_40_markers_M1toM12_OUTPUT_", "",JLidx)
JLidx <- gsub("Extra_40_markers_M13to24_2_", "", JLidx)
mAnnot <- rbind(mAnnot, cbind(JLids, JLidx))
rownames(mAnnot) <- mAnnot[,1]

# Throw away samples without a valid ID
mAnnot <- mAnnot[-which(is.na(mAnnot[,2])),]

# Keep only samples which we have annotated
cvcf <- cvcf[,which(colnames(cvcf) %in% rownames(mAnnot))]
mAnnot <- mAnnot[colnames(cvcf),]
uniqueSamples <- unique(mAnnot[,2])

# Merge all data in one big matrix
monsterplex <- matrix(NA, nrow(cvcf), length(uniqueSamples), dimnames=list(rownames(cvcf), uniqueSamples))

for(x in 1:nrow(mAnnot)){
  cvcfID <- mAnnot[x,1] # ID in cvcf
  umID <- mAnnot[x,2] # UM-Het3 ID in monsterplex
  hasGTdata <- which(!is.na(cvcf[,cvcfID]))
  monsterplex[hasGTdata, umID] <- cvcf[hasGTdata, cvcfID] # Just copy
}

setwd("C:/Github/UM-HET3/files/phenocovar")
covar <- read.csv("ITP_covar_21_July_2020.csv", row.names=1, na.strings=c("U", "x", "", "NA"))
phenos <- read.csv("ITP_phenos_21_July_2020.csv", row.names=1, na.strings=c("U", "x", "", "NA"))

setwd("C:/Github/UM-HET3/files/monsterplex")

monsterplex <- monsterplex[, which(colnames(monsterplex) %in% rownames(covar))]
all(colnames(monsterplex) %in% rownames(phenos))

allP <- cbind(covar[colnames(monsterplex),], phenos[colnames(monsterplex),])
#allP[, "BirthDate_HET3_ITP"] <- as.Date(allP[, "BirthDate_HET3_ITP"]-2, origin = "1900-01-01")
#allP[, "DeathDate_HET3_ITP"] <- as.Date(allP[, "DeathDate_HET3_ITP"]-2, origin = "1900-01-01")

allP <- allP[colnames(monsterplex),c("Sex", "Longevity_HET3_ITP", "Drug_Site", "Site", "Cohort.Year", "BodyWeight_HET3_ITP_6m", "BodyWeight_HET3_ITP_12m", "BodyWeight_HET3_ITP_18m", "BodyWeight_HET3_ITP_24m")]
allP[,"Drug_Site"] <- unlist(lapply(strsplit(allP[,"Drug_Site"], "_"),"[",1))
colnames(allP)[3] <- c("Treatment_Effect")

write.table(monsterplex, "monsterplex.vcfdata.dedup.txt", quote=FALSE, sep="\t")
write.table(map, "monsterplex.snp.annot.txt", quote=FALSE, sep="\t")
write.table(allP, "monsterplex.ind.annot.txt", quote=FALSE, sep="\t")

