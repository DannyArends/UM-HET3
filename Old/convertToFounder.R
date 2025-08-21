
setwd("D:/Edrive/Mouse/ITP Data/VCFS")

f2 <- colnames(cvcf)[7:ncol(cvcf)]

autosomal <- which(cvcf[, "X.CHROM"] != "X")
cgts <- cvcf[autosomal,]
fvcfA <- fvcf[autosomal,]
undecomposable <- c()
for(x in 1:nrow(fvcfA)){
  if(sum(fvcfA[x, 7:10] == "1/1") == 1){ # only one founder has a SNPs
    if(fvcfA[x, "BALB_cJ"] == "1/1"){
      cgts[x, f2[which(cgts[x, f2] == "0/0")]] <-  "B?"
      cgts[x, f2[which(cgts[x, f2] == "0/1")]] <-  "A?"
      cgts[x, f2[which(cgts[x, f2] == "1/1")]] <-  "XX"
    }
    if(fvcfA[x, "C3H_HeJ"] == "1/1"){
      cgts[x, f2[which(cgts[x, f2] == "0/0")]] <-  "?D"
      cgts[x, f2[which(cgts[x, f2] == "0/1")]] <-  "?C"
      cgts[x, f2[which(cgts[x, f2] == "1/1")]] <-  "XX"
    }
    if(fvcfA[x, "DBA_2J"] == "1/1"){
      cgts[x, f2[which(cgts[x, f2] == "0/0")]] <-  "?C"
      cgts[x, f2[which(cgts[x, f2] == "0/1")]] <-  "?D"
      cgts[x, f2[which(cgts[x, f2] == "1/1")]] <-  "XX"
    }
  }
  if(sum(fvcfA[x, 7:10] == "1/1") == 2){ # 2 founders have a SNPs
    if(fvcfA[x, "BALB_cJ"] == "1/1") {
      if(fvcfA[x, "C3H_HeJ"] == "1/1"){
        cgts[x, f2[which(cgts[x, f2] == "0/0")]] <-  "BD"
        cgts[x, f2[which(cgts[x, f2] == "0/1")]] <-  "??"
        cgts[x, f2[which(cgts[x, f2] == "1/1")]] <-  "AC"
      }
      if(fvcfA[x, "DBA_2J"] == "1/1"){
        cgts[x, f2[which(cgts[x, f2] == "0/0")]] <-  "BC"
        cgts[x, f2[which(cgts[x, f2] == "0/1")]] <-  "??"
        cgts[x, f2[which(cgts[x, f2] == "1/1")]] <-  "AD"
      }
    }else{
      cat("Marker", x, "cannot be decomposed into founder genotypes\n")
      undecomposable <- c(undecomposable, x)
    }
  }
  if(sum(fvcfA[x, 7:10] == "1/1") == 3){ # All founders (except for the B6 have a SNP)
    cgts[x, f2[which(cgts[x, f2] == "1/1")]] <-  "A?"
    cgts[x, f2[which(cgts[x, f2] == "0/1")]] <-  "B?"
    cgts[x, f2[which(cgts[x, f2] == "0/0")]] <-  "XX"
  }
}
if(length(undecomposable) > 0) cgts <- cgts[-undecomposable, ]

cat("Autosomal:", nrow(cgts), "undec:", length(undecomposable), "\n")

sexosomal <- which(cvcf[, "X.CHROM"] == "X")
cgtsX <- cvcf[sexosomal,]
fvcfX <- fvcf[sexosomal,]

undecomposable <- c()
for(x in 1:nrow(fvcfX)){
  if(sum(fvcfX[x, 7:10] == "1/1") == 1){ # only one founder has a SNPs
    if(fvcfX[x, "BALB_cJ"] == "1/1"){
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/0")]] <-  "B?"
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/1")]] <-  "A?"
      cgtsX[x, f2[which(cgtsX[x, f2] == "1/1")]] <-  "A-"
    }
    if(fvcfX[x, "C3H_HeJ"] == "1/1"){
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/0")]] <-  "?-"
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/1")]] <-  "?C"
      cgtsX[x, f2[which(cgtsX[x, f2] == "1/1")]] <-  "XX"
    }
    if(fvcfX[x, "DBA_2J"] == "1/1"){
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/0")]] <-  "??"
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/1")]] <-  "XX"
      cgtsX[x, f2[which(cgtsX[x, f2] == "1/1")]] <-  "XX"
    }
  }else if(sum(fvcfX[x, 7:10] == "1/1") == 2){ # Two founders have a SNPs
    if(fvcfX[x, "BALB_cJ"] == "1/1" && fvcfX[x, "C3H_HeJ"] == "1/1"){
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/0")]] <-  "B-"
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/1")]] <-  "??"
      cgtsX[x, f2[which(cgtsX[x, f2] == "1/1")]] <-  "AC"
    }
    if(fvcfX[x, "C57BL_6NJ"] == "1/1" && fvcfX[x, "C3H_HeJ"] == "1/1"){
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/0")]] <-  "XX"
      cgtsX[x, f2[which(cgtsX[x, f2] == "0/1")]] <-  "??"
      cgtsX[x, f2[which(cgtsX[x, f2] == "1/1")]] <-  "BC"
    }
  }else{
    undecomposable <- c(undecomposable, x)
  }
}
if(length(undecomposable) > 0) cgtsX <- cgtsX[-undecomposable, ]

cat("Sexosomal:", nrow(cgtsX), "undec:", length(undecomposable), "\n")

cgts <- rbind(cgts, cgtsX)
cat("Combined:", nrow(cgts), "\n")

