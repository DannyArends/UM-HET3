#
# CTL_Table.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Compute correlation deltas for each SOMA locus as well as effect sizes as: SomaEffects In DaysPGram
#

library(qtl)

source("ActuarialMapping/adjustXprobs.R")
mcross <- read.cross(format="csvr", file="DataSet/um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

#Sample names
snames <- as.character(pull.pheno(mcross)[, "GenoID"])

#42 days
bw <- read.csv("DataSet/bodyweights/ITP_50601.csv", header = TRUE, comment.char = "#", skip=11, row.names=2,na.strings = c("NA", "", "x"))
bw <- bw[which(bw[, "DA2024"] == 1),]

#12 months
trait <- read.csv("DataSet/bodyweights/ITP_10003.csv", header = TRUE, comment.char = "#", skip=10, row.names=2,na.strings = c("NA", "", "x"))
trait <- trait[which(trait[, "DA2024"] == 1),]

adjustPHE <- function(cdata, days = 0, column = "bw42", out = c("is42", "adjBw42", "adjLs42")){
  iz <- cdata[, "longevity"] >= days & !is.na(cdata[, column])
  mdata <- cdata[which(iz),]

  cdata[, out[1]] <<- iz

  lm.null <- lm(as.formula(paste0(column, " ~ sex + site + cohort + treatment")), data = mdata)
  cdata[which(iz), out[2]] <<- round(as.numeric(coef(lm.null)["(Intercept)"]) + residuals(lm.null), 2)

  lm.null.long <- lm(longevity ~ sex + site + cohort + treatment, data = mdata)
  cdata[which(iz), out[3]] <<- round(as.numeric(coef(lm.null.long)["(Intercept)"]) + residuals(lm.null.long),2)
}

gtsp <- pull.genoprob(mcross)
cdata <- data.frame(longevity = as.numeric(pull.pheno(mcross)[, "longevity"]), 
                    sex = as.numeric(pull.pheno(mcross)[, "sex"]), 
                    site = as.factor(pull.pheno(mcross)[, "site"]),
                    bw42 = as.numeric(bw[snames, "Value"]),
                    bw6 = as.numeric(pull.pheno(mcross)[, "bw6"]),
                    bw12 = as.numeric(trait[snames, "Value"]),
                    bw18 = as.numeric(pull.pheno(mcross)[, "bw18"]),
                    bw24 = as.numeric(pull.pheno(mcross)[, "bw24"]),
                    is42 = NA,
                    adjBw42 = NA,
                    adjLs42 = NA, 

                    is6 = NA,
                    adjBw6 = NA,
                    adjLs6 = NA, 

                    is12 = NA,
                    adjBw12 = NA,
                    adjLs12 = NA, 

                    is18 = NA,
                    adjBw18 = NA,
                    adjLs18 = NA, 

                    is24 = NA,
                    adjBw24 = NA,
                    adjLs24 = NA, 
                    cohort = as.factor(pull.pheno(mcross)[, "cohort"]), 
                    treatment = as.factor(pull.pheno(mcross)[, "treatment"]))

adjustPHE(cdata, 0, "bw42", c("is42", "adjBw42", "adjLs42"))
adjustPHE(cdata, 185, "bw6", c("is6", "adjBw6", "adjLs6"))
adjustPHE(cdata, 365, "bw12", c("is12", "adjBw12", "adjLs12"))
adjustPHE(cdata, 550, "bw18", c("is18", "adjBw18", "adjLs18"))
adjustPHE(cdata, 725, "bw24", c("is24", "adjBw24", "adjLs24"))

all <- c("1_3010272","1_86216552","2_13600088", "2_60201233","2_161871392","3_87974845","3_159581164","4_30761996","4_107374161","6_8006720",
         "6_138658041","7_16072018","7_120086292","8_71684276","8_126505019","9_51116640","10_18144599","11_97448477","12_71677220",
         "12_118179607","13_19367506","13_98521647","14_30957748", "14_101437466","15_3288506","16_75758401","17_26542857",
         "18_52488251","19_3403302","19_53851357")
names(all) <- c("Soma1a","Soma1b","Soma2a","Soma2b","Soma2c","Soma3a","Soma3b","Soma4a","Soma4b","Soma6a","Soma6b","Soma7a","Soma7b","Soma8a",
                "Soma8b","Soma9a","Soma10a","Soma11a","Soma12a","Soma12b","Soma13a","Soma13b","Soma14a","Soma14b","Soma15a","Soma16a",
                "Soma17a","Soma18a","Soma19a","Soma19b")

tp <- c(6, 24, 6, 42, 42,18,6,42,42,18,24,12,18,12,42,6,6,12,6,24,42,42,18,42,6,18,24,42,24,24)
sex <- c("C","M","M", "C","M","M","M","M","C","C","F","M","M","M","C","C","C","F","F","M","M","F","M","M","M","C","F","C","M","M")

adata <- cbind(t(t(cbind(names(all), all)[,-1])), tp, sex)

### All SOMA
cat("", file = "DataSet/output/SomaEffectsInDaysPGram_March27.txt", append = FALSE)
for(x in 1:nrow(adata)){
  mp <- gtsp[, grep(adata[x, 1], colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  cat(paste0(adata[x, 1], " "), file = "DataSet/output/SomaEffectsInDaysPGram_March27.txt", append = TRUE)
  for(sex in list(0,1,c(0,1))){

    CH <- which(gts == "AC" & cdata[, "sex"] %in% sex)
    BH <- which(gts == "BC" & cdata[, "sex"] %in% sex)
    CD <- which(gts == "AD" & cdata[, "sex"] %in% sex)
    BD <- which(gts == "BD" & cdata[, "sex"] %in% sex)

    bwC <- paste0("adjBw", adata[x, "tp"])
    lsC <- paste0("adjLs", adata[x, "tp"])
    if(length(CH) > 100) cCH <- lm(cdata[CH, lsC] ~ cdata[CH, bwC]+1);
    if(length(BH) > 100) cBH <- lm(cdata[BH, lsC] ~ cdata[BH, bwC]+1);
    if(length(CD) > 100) cCD <- lm(cdata[CD, lsC] ~ cdata[CD, bwC]+1);
    if(length(BD) > 100) cBD <- lm(cdata[BD, lsC] ~ cdata[BD, bwC]+1);
    cat(round(coef(cCH)[2],1),round(coef(cBH)[2],1),round(coef(cCD)[2],1),round(coef(cBD)[2],1), " ", 
        file = "DataSet/output/SomaEffectsInDaysPGram_March27.txt", append = TRUE)
  }
  cat("\n", file = "DataSet/output/SomaEffectsInDaysPGram_March27.txt", append = TRUE)
}



### All SOMA
cat("", file = "DataSet/output/SomaEffectsCor_March27.txt", append = FALSE)
for(x in 1:nrow(adata)){
  mp <- gtsp[, grep(adata[x, 1], colnames(gtsp))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1,function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))

  cat(paste0(adata[x, 1], " "), file = "DataSet/output/SomaEffectsCor_March27.txt", append = TRUE)
  for(sex in list(0,1,c(0,1))){

    CH <- which(gts == "AC" & cdata[, "sex"] %in% sex)
    BH <- which(gts == "BC" & cdata[, "sex"] %in% sex)
    CD <- which(gts == "AD" & cdata[, "sex"] %in% sex)
    BD <- which(gts == "BD" & cdata[, "sex"] %in% sex)

    bwC <- paste0("adjBw", adata[x, "tp"])
    lsC <- paste0("adjLs", adata[x, "tp"])
    if(length(CH) > 100) cCH <- cor(cdata[CH, lsC], cdata[CH, bwC], method = "spearman", use="pair");
    if(length(BH) > 100) cBH <- cor(cdata[BH, lsC], cdata[BH, bwC], method = "spearman", use="pair");
    if(length(CD) > 100) cCD <- cor(cdata[CD, lsC], cdata[CD, bwC], method = "spearman", use="pair");
    if(length(BD) > 100) cBD <- cor(cdata[BD, lsC], cdata[BD, bwC], method = "spearman", use="pair");
    cat(round(cCH,2),round(cBH,2),round(cCD,2),round(cBD,2), " ", file = "DataSet/output/SomaEffectsCor_March27.txt", append = TRUE)
  }
  cat("\n", file = "DataSet/output/SomaEffectsCor_March27.txt", append = TRUE)
}

