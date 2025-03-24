mm <- c("6_54992703", "6_62673982", "6_81943533", "6_83032867", "6_91043806", "6_92807055", "6_97112728", "6_102489070", "6_108075853", "6_118542923", "6_119943820", "6_122290952")

mm <- c("8_55830812", "8_65000838", "8_70484370", "8_72345998", "8_81997313", "8_86034949", "8_91341862", "8_95039510", "8_107119885", "8_111319339", "8_124035048", "8_124525717")

mm <- c("2_60201233", "2_62606286", "2_65326540", "2_71846197", "8_81997313", "2_83915232", "2_92975818", "2_99777113", "2_105696373", "2_110880795", "2_122685194", "2_125083168")


setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)
markers <- colnames(pull.geno(mcross))

mm <- c("6_81943533", "6_83032867", "6_91043806", "6_92807055", "6_97112728", "6_102489070", "6_108075853")

GTS <- c()
for(marker in mm){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  ## Maternal / Paternal
  #gts[gts == "AC"] <- "A"; gts[gts == "AD"] <- "A"
  #gts[gts == "BC"] <- "B"; gts[gts == "BD"] <- "B"

  ## Maternal / Paternal
  gts[gts == "AC"] <- "C"; gts[gts == "AD"] <- "D"
  gts[gts == "BC"] <- "C"; gts[gts == "BD"] <- "D"


  GTS <- cbind(GTS, gts)
}

table(unlist(lapply(apply(GTS, 1, table), length)))

plot(c(1, length(markers)), y = c(0, 0.4), t = "n")
nn <- c()
for(x in 1:(length(markers) - 8)){
  GTS <- c()
  for(marker in markers[x:(x+8)]){
    mp <- gtsp[, grep(marker, colnames(gtsp))]
    gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
      if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
    }))
    ## Maternal
    #gts[gts == "AC"] <- "A"; gts[gts == "AD"] <- "A"
    #gts[gts == "BC"] <- "B"; gts[gts == "BD"] <- "B"

    ## Paternal
    gts[gts == "AC"] <- "C"; gts[gts == "AD"] <- "D"
    gts[gts == "BC"] <- "C"; gts[gts == "BD"] <- "D"

    ## BXD
    #gts[gts == "AC"] <- "C"; gts[gts == "AD"] <- "D"
    #gts[gts == "BC"] <- "C"; gts[gts == "BD"] <- "D"

    GTS <- cbind(GTS, gts)
  }
  tbl <- table(unlist(lapply(apply(GTS, 1, table), length)))
  n <- tbl["2"] / (tbl["1"]+tbl["2"])
  points(x, n, col = "green")
  nn <- c(nn, n)
}


cbind(markers[as.numeric(which(nn < 0.001))], markers[8+as.numeric(which(nn < 0.001))])

aa <- scanone(mcross)
p <- c()
for(x in 1:(nrow(aa)-8)){
  p <- c(p, aa[x+8,2]-aa[x,2])
}

555
