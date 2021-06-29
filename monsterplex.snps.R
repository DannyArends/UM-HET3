#Create the children (cvcf file)

setwd("D:/Edrive/Mouse/ITP Data/VCFSmay2021")

cvcf <- c()
for(chr in c(1:19, "X")){
  mvcf <- read.table(paste0("snps-filtered_chr", chr,".vcf"), sep="\t", comment.char="", skip=95, header=TRUE, colClasses="character")
  complexSNP <- grep(",", mvcf[,5])
  if(length(complexSNP) > 0) mvcf <- mvcf[-complexSNP,] # Take only simple SNPs (2 alleles)
  mvcf <- mvcf[which(nchar(mvcf[,4]) == 1 & nchar(mvcf[,5]) == 1),] # Take only SNPs

  for (x in 1:nrow(mvcf)) {
    mvcf[x, 10:ncol(mvcf)] <- unlist(lapply(strsplit(as.character(mvcf[x,10:ncol(mvcf)]), ":"),"[",1))
    nas <- which(mvcf[x, 10:ncol(mvcf)] == "./.")
    mvcf[x, (10:ncol(mvcf))[nas]] <- NA
  }
  mvcf <- mvcf[,-c(3,7,8,9)]
  cvcf <- rbind(cvcf, mvcf)
  cat("Loaded chr", chr, "\n")
}

cvcf[,1] <- gsub("chr", "", cvcf[,1])
write.table(cvcf, file="cvcfAll.txt",sep="\t", quote=FALSE, row.names=FALSE)

write.table(cvcf[,1:2], file = "regionsITP.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

setwd("D:/Edrive/Mouse/ITP Data/VCFSmay2021")
fvcfRaw <- read.table("Founders.snps-filtered.vcf", sep="\t", comment.char="", skip=95, header=TRUE)

complexSNP <- grep(",", fvcfRaw[,5])
if(length(complexSNP) > 0) fvcfRaw <- fvcfRaw[-complexSNP,] # Take only simple SNPs (2 alleles)
fvcfRaw <- fvcfRaw[which(nchar(fvcfRaw[,4]) == 1 & nchar(fvcfRaw[,5]) == 1),] # Take only SNPs

for (x in 1:nrow(fvcfRaw)) {
  fvcfRaw[x, 10:ncol(fvcfRaw)] <- unlist(lapply(strsplit(as.character(fvcfRaw[x,10:ncol(fvcfRaw)]), ":"),"[",1))
  nas <- which(fvcfRaw[x, 10:ncol(fvcfRaw)] == "./.")
  fvcfRaw[x, (10:ncol(fvcfRaw))[nas]] <- NA
}
fvcfRaw <- fvcfRaw[,-c(3,7,8,9)]
colnames(fvcfRaw)[6] <- "BALBBYJ"
write.table(fvcfRaw, file="fvcfAll.txt",sep="\t", quote=FALSE, row.names=FALSE)
