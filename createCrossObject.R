setwd("C:/Github/UM-HET3/files")

map <- read.table("merged/map.gts4way.Juli2021.txt", sep="\t")
ind <- read.table("phenotypes/HET3-ITP_traitsGN2_Sep2021.csv", sep=",", skip=11, header=TRUE,na.strings="x")
ind <- ind[-c(1:20),]
ind <- ind[-grep("_SE", ind[,1]),]
rownames(ind) <- ind[,1]
ind <- ind[,-1]

# Animals marked as REMOVED, but genotyped
removed <- c("UT04802", "UT04803", "UT04501", "UT04502", "UT04539", "UT04540", 
             "UT04541", "UT04542", "UT04725", "UT04726", "UT04727", "UT04728",
             "UM39151", "UM40061", "UM40467", "UM41174", "UM41419", "UM41426", "UM41993", "UM42007", "UM43218")

GTS <- read.table("ITP_6480x893_MatPat_Sep21.txt", sep="\t")
mGTS <- read.table("ITP_6480x486_Maternal_Sep21.txt", sep="\t")
pGTS <- read.table("ITP_6480x396_Paternal_Sep21.txt", sep="\t")

individuals <- colnames(GTS)[-c(1:2)]
ind <- ind[individuals,]
ind <- cbind(ind, sex = 1 - as.numeric(ind[,"Sex"]))
ind <- ind[-which(ind[, "Tx_Group"] == "2.0"),]
ind <- ind[-which(rownames(ind) %in% removed),]

dim(ind)
GTS <- GTS[, c("Chr","Position", rownames(ind))]
mGTS <- mGTS[, c("Chr","Position", rownames(ind))]
pGTS <- pGTS[, c("Chr","Position", rownames(ind))]

write.table(GTS, file="ITP_6438x893_MatPat_Sep21.txt", sep = "\t", quote=FALSE,na="")
write.table(mGTS, file="ITP_6438x486_Maternal_Sep21.txt", sep = "\t", quote=FALSE,na="")
write.table(pGTS, file="ITP_6438x396_Paternal_Sep21.txt", sep = "\t", quote=FALSE,na="")

indM <- cbind(sex = ind[, "sex"], 
              treatment = ind[, "Tx_Group"], 
              site = ind[, "Site"], 
              cohort = ind[, "Year"], 
              bw6 = ind[, "ITP_Weight6m"], bw12 = ind[, "ITP_Weight12m"], bw18 = ind[, "Body_18m"], bw24 = ind[, "BodyWeight_HET3_ITP_24m"],
              longevity = ind[, "Longevity_HET3_ITP"]
           )

# Write out our cross object for QTL mapping
write.table(cbind(NA,NA, t(cbind(Individual = rownames(indM), indM))), file = "um-het3-rqtl.csvr", col.names = FALSE, sep = ",", quote=FALSE, na = "")
write.table(rbind(GenoID = c(NA, NA, colnames(GTS)[-c(1,2)]), 
                  cbind(Chr = GTS[,"Chr"], Mb = as.numeric(GTS[,"Position"]) / 1000000, GTS[, -c(1,2)])
                 ), file = "um-het3-rqtl.csvr", col.names = FALSE, sep = ",", append=TRUE, quote=FALSE, na="")

                 
## find gaps larger than 10 mb

GTS[sort(c(which(diff(GTS[,2]) > 10000000), which(diff(GTS[,2]) > 10000000)+1)),1:2]
mGTS[sort(c(which(diff(mGTS[,2]) > 15000000), which(diff(mGTS[,2]) > 15000000)+1)),1:2]
pGTS[sort(c(which(diff(pGTS[,2]) > 15000000), which(diff(pGTS[,2]) > 15000000)+1)),1:2]