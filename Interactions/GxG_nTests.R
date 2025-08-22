#
# GxGntests.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Code to compute the exact number of tests performed when doing the G x G scan across all time points
#

library(RColorBrewer)
library(svglite)
library(vioplot)

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")


allV <- c("1_3010274", "1_24042124", "1_121483290", "1_167148678", "2_89156987", "2_112255823", "2_148442635", "3_83354281", 
         "4_52524395", "4_154254581", "5_67573068", "6_93680853", "6_132762500", "9_34932404", "9_104091597", "9_124056586", 
         "10_72780332", "11_6599922", "11_82176894", "11_113729074", "12_112855820", "13_83858506", "14_78415875", "14_101437466", 
         "15_74248242", "15_99306167", "17_32883804", "18_52488251", "X_36008085")

names(allV) <- c("Vita1a", "Vita1b", "Vita1c", "Vita1d", "Vita2a", "Vita2b", "Vita2c", "Vita3a", "Vita4a", "Vita4b", "Vita5a", "Vita6a", "Vita6b", 
                "Vita9a", "Vita9b", "Vita9c", "Vita10a", "Vita11a", "Vita11b", "Vita11c", "Vita12a", "Vita13a", "Vita14a", "Vita14b", "Vita15a", 
                "Vita15b", "Vita17a", "Vita18a", "VitaXa")

allS <- c("1_3010274","1_86216552","2_13600088", "2_60201233","2_161871392","3_87974845","3_159581164","4_30761996","4_107374161","6_8006720",
         "6_138658041","7_16072018","7_120086292","8_71684276","8_126505019","9_51116640","10_18144599","11_97448477","12_71677220",
         "12_118179607","13_19367506","13_98521647","14_30957748", "14_101437466","15_3288506","16_75758401","17_26542857",
         "18_52488251","19_3403302","19_53851357")
names(allS) <- c("Soma1a","Soma1b","Soma2a","Soma2b","Soma2c","Soma3a","Soma3b","Soma4a","Soma4b","Soma6a","Soma6b","Soma7a","Soma7b","Soma8a",
                "Soma8b", "Soma9a", "Soma10a","Soma11a","Soma12a","Soma12b","Soma13a","Soma13b","Soma14a","Soma14b","Soma15a","Soma16a",
                "Soma17a","Soma18a","Soma19a","Soma19b")

all <- c(allV, allS)
all <- all[-which(duplicated(all))]

pos <- as.numeric(unlist(lapply(strsplit(all, "_"),"[", 2)))
pos.ix <- sort(pos, index.return = TRUE)$ix

all <- all[pos.ix]
chrs <- unlist(lapply(strsplit(all, "_"),"[",1))
chrs[chrs == "X"] <- "20"
chrs.ix <- sort(as.numeric(chrs), index.return = TRUE)$ix

all <- all[chrs.ix]

int.f <- c()
int.m <- c()
#for(tp in c("42", "365", "740", "905")){

tp <- 42

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny/x_GxG_Vita_Soma_See_Extended_Data_6")
lodM.m <- read.table(paste0("vita_soma_interactions_2way_males_tp",tp,".txt"), sep = "\t")
for(x in 1:ncol(lodM.m)){
for(y in 1:x){
    lodM.m[y,x] <- lodM.m[x,y]
  }
} 
lodM.m <- lodM.m[names(all), names(all)]
lodM.f <- read.table(paste0("vita_soma_interactions_2way_females_tp",tp,".txt"), sep = "\t")
for(x in 1:ncol(lodM.f)){
for(y in 1:x){
    lodM.f[y,x] <- lodM.f[x,y]
  }
} 
lodM.f <- lodM.f[names(all), names(all)]

colz.c <- colorRampPalette(c("white", "lightskyblue3"))(6)
colz.c <- c(colz.c, rep("lightskyblue3", 40))

colz.c2 <- colorRampPalette(c("white", "plum2"))(6)
colz.c2 <- c(colz.c2, rep("plum2", 40))

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny/")

# TODO: Make these images SEX-Specific

#pdf(paste0(paste0("GxG_tp",tp,"_Vita_Soma_Big.pdf")), width = 36, height = 36)
#op <- par(mar = c(5,5,5,5))
#plot(c(0.5, length(all)+0.5), c(0.5, length(all)+0.5), t = "n", xaxt='n', yaxt='n', xlab="",ylab="", xaxs="i", yaxs="i",bty="n")
#axis(1, at = 1:length(all), rownames(lodM.m[names(all),])[1:length(all)],las=2)
#axis(2, at = 1:length(all), rev(rownames(lodM.m[names(all),])[1:length(all)]),las=2)
#axis(3, at = 1:length(all), rownames(lodM.m[names(all),])[1:length(all)],las=2)
#axis(4, at = 1:length(all), rev(rownames(lodM.m[names(all),])[1:length(all)]),las=2)
nTests <- 0
nSig.m <- c(0, 0, 0)
nSig.f <- c(0, 0, 0)

allLods.m <- c()
allLods.f <- c()

for(x in 1:length(all)){
  for(y in 1:length(all)){
    #cat("x:",x,"y:",y, "\n")
    xp <- x
    name1 <- names(all)[x]
    yp <- (length(all)+1) - y
    name2 <- names(all)[y]
    m1 <- gsub("Soma", "", gsub("Vita", "", name1))
    m2 <- gsub("Soma", "", gsub("Vita", "", name2))
    m1 <- substr(m1, 1, nchar(m1)-1)
    m2 <- substr(m2, 1, nchar(m2)-1)
    if(m1 != m2){
      if(x > y && !is.na(lodM.m[name2,name1])){
        bcol <- "gray"
        if(lodM.m[name2, name1] > 3.8){
          bcol = "black"
          t = 0
          if(grepl("Vita", name2) && grepl("Vita", name1)){ nSig.m[1] <- nSig.m[1] + 1; t = 1; }
          if(grepl("Soma", name2) && grepl("Soma", name1)){ nSig.m[2] <- nSig.m[2] + 1; t = 2; }
          if(grepl("Vita", name2) && grepl("Soma", name1)){ nSig.m[3] <- nSig.m[3] + 1; t = 3; }
          if(grepl("Soma", name2) && grepl("Vita", name1)){ nSig.m[3] <- nSig.m[3] + 1; t = 3; }
          int.m <- rbind(int.m, c(name2, name1, tp, t, round(lodM.m[name2,name1],2)))
        }
        #rect(xp-0.5,yp-0.5, xp+0.5, yp + 0.5, col = colz.c[1+ 0.5 * round(lodM.m[name2,name1])], border = "gray")
        #text(xp, yp, paste0(formatC(lodM.m[name2,name1], digits = 1, format = "f"), ""), col = bcol)
        nTests <- nTests + 1
        allLods.m <- c(allLods.m, lodM.m[name2, name1])
      }
      if(x < y &&  !is.na(lodM.f[name2,name1])){
        bcol <- "gray"
        if(lodM.f[name2,name1] > 3.8){
          bcol = "black"
          if(grepl("Vita", name2) && grepl("Vita", name1)){ nSig.f[1] <- nSig.f[1] + 1; t = 1; }
          if(grepl("Soma", name2) && grepl("Soma", name1)){ nSig.f[2] <- nSig.f[2] + 1; t = 2; }
          if(grepl("Vita", name2) && grepl("Soma", name1)){ nSig.f[3] <- nSig.f[3] + 1; t = 3; }
          if(grepl("Soma", name2) && grepl("Vita", name1)){ nSig.f[3] <- nSig.f[3] + 1; t = 3; }
          int.f <- rbind(int.f, c(name2, name1, tp, t, round(lodM.f[name2, name1],2)))
        }
        #rect(xp-0.5,yp-0.5, xp+0.5, yp + 0.5, col = colz.c2[1+ 0.5 * round(lodM.f[name2,name1])], border = "gray")
        #text(xp, yp, paste0(formatC(lodM.f[name2,name1], digits = 1, format = "f"), ""), col = bcol)
        allLods.f <- c(allLods.f, lodM.f[name2, name1])
      }
    }
  }
}

