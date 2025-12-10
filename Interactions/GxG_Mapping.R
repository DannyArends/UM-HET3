#
# GxG_Mapping.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Mapping GxG interactions between significant Vita and Soma loci at 4 different timepoints (42, 365, 740, and 905 days)
#

allV <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085")

names(allV) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa")

allS <- c("1_3010272","1_86216552","2_13600088", "2_60201233","2_161871392","3_87974845","3_159581164","4_30761996","4_107374161","6_8006720",
         "6_138658041","7_16072018","7_120086292","8_71684276","8_126505019","9_51116640","10_18144599","11_97448477","12_71677220",
         "12_118179607","13_19367506","13_98521647","14_30957748", "14_101437466","15_3288506","16_75758401","17_26542857",
         "18_52488251","19_3403302","19_53851357")
names(allS) <- c("Soma1a","Soma1b","Soma2a","Soma2b","Soma2c","Soma3a","Soma3b","Soma4a","Soma4b","Soma6a","Soma6b","Soma7a","Soma7b","Soma8a",
                "Soma8b", "Soma9a", "Soma10a","Soma11a","Soma12a","Soma12b","Soma13a","Soma13b","Soma14a","Soma14b","Soma15a","Soma16a",
                "Soma17a","Soma18a","Soma19a","Soma19b")

all <- c(allV, allS)

library(qtl)

source("ActuarialMapping/adjustXprobs.R")

# Read cross object
mcross <- read.cross(format="csvr", file="DataSet/um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)

sexes <- c(0, 1)
names(sexes) <- c("females", "males")
timepoints <- c(42, 365, 740, 905)

## Combined
for(tp in timepoints){
  gtsp <- pull.genoprob(mcross)
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "longevity"] >= tp & cdata[, "sex"] == "0")

  cdata <- cdata[idx, ]
  gtsp <- gtsp[idx, ]


  # Compute GxG LOD scores for males and females combined

  lodM <- matrix(NA, length(all), length(all), dimnames = list(names(all), names(all)))
  for(m1 in 1:length(all)){
    lodR <- c()
    for(m2 in m1:length(all)){
      mp1 <- gtsp[, grep(all[m1], colnames(gtsp))]
      mp2 <- gtsp[, grep(all[m2], colnames(gtsp))]
      lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = cdata)
      lm.alt <- lm(longevity ~  site + cohort + treatment + mp1 + mp2 + 0, data = cdata)
      lm.alt2 <- lm(longevity ~  site + cohort + treatment + mp1 * mp2 + 0, data = cdata)
      n <- sum(!is.na(lm.alt$resid))
      lodM[names(all)[m1], names(all)[m2]] <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
      lodM[names(all)[m2], names(all)[m1]] <- (n/2) * log10(sum(lm.alt$resid^2) / sum(lm.alt2$resid^2))
      if(m1 == m2) lodM[m1,m2] <- NA
    }
  }

  #cat(names(coef(lm.alt2)), "\n", coef(lm.alt2), "\n")
  rownames(lodM) <- names(all)
  colnames(lodM) <- names(all)

  library(RColorBrewer)
  colz <- c("white", brewer.pal(4,"PuRd"))

  image(1:length(all), 1:length(all), lodM, xaxt="n", yaxt="n", xlab="", ylab="", breaks = c(0, 3, 6, 9, 12, 100), col=colz)
  axis(1, at = 1:length(all), names(all), las=2)
  axis(2, at = 1:length(all), names(all), las=2)

  write.table(lodM, paste0("DataSet/output/vita_soma_interactions_2way_combined_tp", tp,".txt"), sep = "\t", quote=FALSE)
}

### Males and Females
for(sex in names(sexes)){
  for(tp in timepoints){
    gtsp <- pull.genoprob(mcross)
    cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                        sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                        site = as.factor(pull.pheno(mcross)[, "site"]),
                        cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                        treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

    idx <- which(cdata[, "sex"] == sexes[sex] & cdata[, "longevity"] >= tp)

    cdata <- cdata[idx, ]
    gtsp <- gtsp[idx, ]


    # Compute GxG LOD scores for males or females

    lodM <- matrix(NA, length(all), length(all), dimnames = list(names(all), names(all)))
    for(m1 in 1:length(all)){
      lodR <- c()
      for(m2 in m1:length(all)){
        mp1 <- gtsp[, grep(all[m1], colnames(gtsp))]
        mp2 <- gtsp[, grep(all[m2], colnames(gtsp))]
        lm.null <- lm(longevity ~ site + cohort + treatment + 0, data = cdata)
        lm.alt <- lm(longevity ~  site + cohort + treatment + mp1 + mp2 + 0, data = cdata)
        lm.alt2 <- lm(longevity ~  site + cohort + treatment + mp1 * mp2 + 0, data = cdata)
        n <- sum(!is.na(lm.alt$resid))
        lodM[names(all)[m1], names(all)[m2]] <- (n/2) * log10(sum(lm.null$resid^2) / sum(lm.alt$resid^2))
        lodM[names(all)[m2], names(all)[m1]] <- (n/2) * log10(sum(lm.alt$resid^2) / sum(lm.alt2$resid^2))
        if(m1 == m2) lodM[m1,m2] <- NA
      }
    }
    rownames(lodM) <- names(all)
    colnames(lodM) <- names(all)

    library(RColorBrewer)
    colz <- c("white", brewer.pal(4,"PuRd"))

    image(1:length(all), 1:length(all), lodM, xaxt="n", yaxt="n", xlab="", ylab="", breaks = c(0, 3, 6, 9, 12, 100), col=colz)
    axis(1, at = 1:length(all), names(all), las=2)
    axis(2, at = 1:length(all), names(all), las=2)

    write.table(lodM, paste0("DataSet/output/vita_soma_interactions_2way_",sex,"_tp", tp,".txt"), sep = "\t", quote=FALSE)
  }
}

### DO the green dots
combos <- rbind(
            c("Vita2a", "Vita17a"),
            c("Vita11a", "Vita14a"),
            c("Vita4a", "Vita10a"),
            c("Vita10a", "Vita11a"),
            c("Vita3a", "Vita6b"),
            c("Vita11b", "Vita13a"),
            c("Vita11a", "Vita13a"),
            c("Vita11a", "Vita15a"),
            c("Vita2c", "Vita18a"),
            c("Vita3b", "Vita4a"),
            c("Vita1b", "Vita10a"),
            c("Vita1c", "Vita13a"),
            c("Vita6a", "Vita12a"))

for(x in 1:nrow(combos)){
  m1 <- combos[x,1]
  m2 <- combos[x,2]
  tp <- timepoints[4]

  gtsp <- pull.genoprob(mcross)
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "sex"] == 0 & cdata[, "longevity"] > tp)

  cdata <- cdata[idx, ]
  gtsp <- gtsp[idx, ]
  mp1 <- gtsp[, grep(all[m1], colnames(gtsp))]
  mp2 <- gtsp[, grep(all[m2], colnames(gtsp))]
  lm.F <- lm(longevity ~  site + cohort + treatment + mp1[,-4] * mp2[,-4] + 0, data = cdata)
  in.F <- coef(lm.F)[19:27]


  gtsp <- pull.genoprob(mcross)
  cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                      sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                      site = as.factor(pull.pheno(mcross)[, "site"]),
                      cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                      treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

  idx <- which(cdata[, "sex"] == 1 & cdata[, "longevity"] > tp)

  cdata <- cdata[idx, ]
  gtsp <- gtsp[idx, ]
  mp1 <- gtsp[, grep(all[m1], colnames(gtsp))]
  mp2 <- gtsp[, grep(all[m2], colnames(gtsp))]
  lm.M <- lm(longevity ~  site + cohort + treatment + mp1[,-4] * mp2[,-4] + 0, data = cdata)
  in.M <- coef(lm.M)[19:27]

  ii <- which(abs(in.F) > 5 & abs(in.M) > 5)
  aa <- cbind(in.F[ii], in.M[ii], sign(in.F[ii]) == sign(in.M[ii]))

  cat(m1,m2, round(sum(aa[,3]) / length(ii),2), "\n")
}

