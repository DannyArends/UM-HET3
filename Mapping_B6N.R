setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/")

markers <- read.table("merged/all.vcf.sorted.txt", sep = "\t")
map <- read.table("genetic_map.txt", sep = "\t")
pheno <- read.table("merged/ind.gts4way.txt", sep = "\t")
parental <- read.table("merged/fvcfAll.txt", sep = "\t")

dim(markers)

markers <- markers[rownames(map), ]
parental <- parental[rownames(map), ]

pheno <- pheno[which(rownames(pheno) %in% colnames(markers)),]

parental <- parental[which(apply(apply(parental,1,is.na),2,sum) == 0),]
markers <- markers[rownames(parental), rownames(pheno)]
map <- map[rownames(parental), ]

dim(map)
dim(markers)
dim(parental)
dim(pheno)

mB6 <- rownames(parental)[which(parental[,"BALB_cJ"] == "1/1" & parental[,"C3H_HeJ"] == "1/1" & parental[,"DBA_2J"] == "1/1")]

# Suitable markers
map <- map[mB6,]
markers <- markers[mB6,]

dim(map)
dim(markers)
dim(pheno)

sInd <- names(which(apply(apply(markers,2,is.na),2,sum) != 371))

markers <- markers[,sInd]
pheno <- pheno[sInd,]

msequence <- seq(365, 1100, 15)
lods.cM <- c()
for(x in msequence){
  lods.c <- c()
  for(marker in rownames(markers)){
    cdata <- data.frame(longevity = as.numeric(pheno[, "Longevity_HET3_ITP"]), 
                        sex = as.factor(pheno[, "Sex"]), 
                        site = as.factor(pheno[, "Site"]),
                        cohort = as.factor(pheno[, "Cohort.Year"]), 
                        treatment = as.factor(pheno[, "Treatment_Effect"] != "Cont"),
                        gt = as.factor(as.character(markers[marker, rownames(pheno)]))
                        )
    noMissing <- apply(apply(cdata,1,is.na),2,sum) == 0
    cdata <- cdata[which(cdata[,"longevity"] > x & noMissing),]
    tryCatch({
      lm.null <- lm(longevity ~ sex + site + cohort + treatment + 0, data = cdata)
      lm.alt <- lm(longevity ~ sex + site + cohort + treatment + gt + 0, data = cdata)
      n <- sum(!is.na(lm.alt$resid))
      lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
      lods.c <- c(lods.c, lod)
    }, error = function(x){
      lods.c <<- c(lods.c, NA)
    })

  }
  lods.cM <- rbind(lods.cM, lods.c)
  cat("Done", x, "\n")
}

lods.cM[is.infinite(lods.cM)] <- NA

mmx <- which(apply(t(lods.cM) > 4,1,sum,na.rm=TRUE) > 0)
mL <- apply(lods.cM[,mmx],2,max)

cbind(map[mmx,], mL)



