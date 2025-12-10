#
# stratifiedPerms.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Use simulations to get a threshold for LOD scores, using permutations stratified by sex and year
#

set.seed(1) # Make sure the permutations are re-do-a-ble
for(p in 1:1000){
  pdata <- cdata
  s1 <- which(pdata[,"sex"]==0 & pdata[,"site"]==1)
  s2 <- which(pdata[,"sex"]==0 & pdata[,"site"]==2)
  s3 <- which(pdata[,"sex"]==0 & pdata[,"site"]==3)
  s4 <- which(pdata[,"sex"]==1 & pdata[,"site"]==1)
  s5 <- which(pdata[,"sex"]==1 & pdata[,"site"]==2)
  s6 <- which(pdata[,"sex"]==1 & pdata[,"site"]==3)

  pdata[s1, "longevity"] <- pdata[sample(s1), "longevity"]
  pdata[s2, "longevity"] <- pdata[sample(s2), "longevity"]
  pdata[s3, "longevity"] <- pdata[sample(s3), "longevity"]
  pdata[s4, "longevity"] <- pdata[sample(s4), "longevity"]
  pdata[s5, "longevity"] <- pdata[sample(s5), "longevity"]
  pdata[s6, "longevity"] <- pdata[sample(s6), "longevity"]

  lods.c <- c()
  lm.null <- lm(longevity ~ sex + site + cohort + treatment + 0, data = pdata)
  for(marker in colnames(pull.geno(mcross))){
    mp <- gtsp[, grep(marker, colnames(gtsp))]
    lm.alt <- lm(longevity ~ sex + site + cohort + treatment + mp + 0, data = pdata)
    n <- sum(!is.na(lm.alt$resid))
    lod <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
    lods.c <- c(lods.c, lod)
  }
  p.lods <- c(p.lods, max(lods.c))
  cat("Done p=",p,", lod=", max(lods.c), "\n")
}

