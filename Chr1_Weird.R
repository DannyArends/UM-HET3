#
# Chr1_Weird.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Code to investigating some weird recombination patterns seen on chromosome 1 between 31.062.759 and 34.173.196 bp
#

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)


all <- c("1_31062759","1_34173196")

names(all) <- c("prox", "dist")

genotypes <- c()
for(marker in all){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  genotypes <- cbind(genotypes, gts)
}
colnames(genotypes) <- all

rec <- sum(unlist(lapply(apply(genotypes,1,table),length)) == 2)
nonrec <- sum(unlist(lapply(apply(genotypes,1,table),length)) == 1)
rec / (nonrec+rec)

all <- c("1_34173196", "1_37045066")

names(all) <- c("prox", "dist")

genotypes <- c()
for(marker in all){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  genotypes <- cbind(genotypes, gts)
}
colnames(genotypes) <- all
cor(apply(genotypes,2,function(x){as.numeric(as.factor(x))}), use = "pair")

rec <- sum(unlist(lapply(apply(genotypes,1,table),length)) == 2)
nonrec <- sum(unlist(lapply(apply(genotypes,1,table),length)) == 1)
rec/ (nonrec+rec)


all <- c("1_26682184", "1_31062759")

names(all) <- c("prox", "dist")

genotypes <- c()
for(marker in all){
  mp <- gtsp[, grep(marker, colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  genotypes <- cbind(genotypes, gts)
}
colnames(genotypes) <- all

rec <- sum(unlist(lapply(apply(genotypes,1,table),length)) == 2)
nonrec <- sum(unlist(lapply(apply(genotypes,1,table),length)) == 1)
rec/ (nonrec+rec)

