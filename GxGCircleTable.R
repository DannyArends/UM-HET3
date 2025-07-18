###
### Candidate gene lists, features in Supplemental 2
###

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

gtsp <- pull.genoprob(mcross)

setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
lods.all <- read.table("progressiveMapping_all20D.txt", sep = "\t", check.names = FALSE)
lods.f <- read.table("progressiveMapping_females20D.txt", sep = "\t", check.names = FALSE)
lods.m <- read.table("progressiveMapping_males20D.txt", sep = "\t", check.names = FALSE)

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

mbigm <- c()

for(tp in c("42", "365", "740", "905")){
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

  setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
  map <- read.table("genetic_map.txt", sep = "\t")
  map <- cbind(map, cloc = NA)
  map <- map[-c(893,892),]


  chrs <- c(1:19, "X")
  gap <- 40*1e6
  chrs.length <- as.numeric(unlist(lapply(chrs, function(x){ max(map[map[, "Chr"] == x,"Pos"]) })))
  names(chrs.length) <- chrs
  l.x <- sum(chrs.length) + gap * (length(chrs)-1)


  circlelocations <- function(nt) {
    medpoints <- matrix(nrow = nt, ncol = 2)
    phi <- seq(0, 2 * pi, length = (nt + 1))
    complex.circle <- complex(modulus = 1, argument = phi)
    for (j in 1:nt) {
      medpoints[j, ] <- c(Im(complex.circle[j]), Re(complex.circle[j]))
    }
    medpoints
  }

  drawspline <- function (cn1, cn2, lwd = 1,col="blue",...){
    x <- cbind(cn1[1],0,cn2[1])
    y <- cbind(cn1[2],0,cn2[2])
    r <- xspline(x, y, lty=1, shape=1, lwd=lwd, border=col,...)
  }

  gsize <- (gap * (length(chrs))) / 1e6

  locs <- circlelocations(sum(ceiling(chrs.length / 1e6)) + gsize)
  locs <- cbind(locs, col = "white")

  chr.starts <- c((gap/1e6)/2)
  cp <-  (gap/1e6) /2
  i <- 1
  col <- "black"
  for(l in chrs.length ) {
    p <- ceiling(l / 1e6)
    #cat(cp, " ", p, " ", cp+p, "\n")
    locs[cp:(cp+p), 3] <- col
    cp <- cp + p + (gap/1e6)
    chr.starts <- c(chr.starts, cp)
    i <- i + 1
    if(col == "black"){ 
      col = "orange"
    }else{ 
      col = "black" 
    }
  }
  locs <- data.frame(locs)
  locs[,1] <- as.numeric(locs[,1])
  locs[,2] <- as.numeric(locs[,2])

  names(chr.starts) <- chrs

  chr.mids <- c()
  for(x in 1:(length(chr.starts)-1)){
    chr.mids <- c(chr.mids, chr.starts[x] + ((chr.starts[x+1] - chr.starts[x] - (gap/1e6)) / 2))
    cat(chr.starts[x], " ", chr.starts[x+1], " = ", chr.mids[x], "\n")
  }
  names(chr.mids) <- names(chr.starts)[1:length(chr.mids)]

  toPos <- function(chr, pos){ return(chr.starts[chr] + round(pos / 1e6)) }

  for(x in 1:nrow(map)){
    map[x, "cloc"] = toPos(map[x, "Chr"], map[x, "Pos"])
  }



### MALES
  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/__Arends_Nature_Prep_All_Key_Files/11_FiguresDanny")
  #pdf(paste0("Circle_v2_",tp,".pdf"), width = 20, height = 10)
  #op <- par(mfrow=c(1,2))
  #par(cex=1)
  #par(pch=3)


#plot(locs[,1:2], pch = 19, cex=0.1, col = locs[, "col"],xlab="", ylab="", xaxt="n", yaxt="n", bty="n", main = paste0("Male", tp))
  mcnt <- 0
#text(1.05 * locs[chr.mids, c(1:2)], names(chr.mids))
  for(x in 1:(nrow(lodM.m)-1)){
    l1 <- map[all[rownames(lodM.m)[x]], "cloc"]
    pch <- 1; col = "black";
    if(grepl("Soma", rownames(lodM.m)[x])){ pch <- 20; col = "green";}

    for(y in (x+1):nrow(lodM.m)){
      m1 <- gsub("Soma", "", gsub("Vita", "", rownames(lodM.m)[x]))
      m2 <- gsub("Soma", "", gsub("Vita", "", rownames(lodM.m)[y]))
      m1 <- substr(m1, 1, nchar(m1)-1)
      m2 <- substr(m2, 1, nchar(m2)-1)

      type <- 1
      if(grepl("Soma", rownames(lodM.m)[x]) && grepl("Vita", rownames(lodM.m)[y])) type <- 2
      if(grepl("Vita", rownames(lodM.m)[x]) && grepl("Soma", rownames(lodM.m)[y])) type <- 2
      if(grepl("Soma", rownames(lodM.m)[x]) && grepl("Soma", rownames(lodM.m)[y])) type <- 3

      if(m1 != m2 && lodM.m[x,y] > 2.0){
        cat(all[rownames(lodM.m)[x]], " ",all[rownames(lodM.m)[y]], "\n")
        s <- 1
        if(lodM.m[x,y] > 3.8 && lodM.m[x,y] <= 4.2) s <- 2
        if(lodM.m[x,y] > 4.2 && lodM.m[x,y] <= 4.5) s <- 3
        if(lodM.m[x,y] > 4.5 ) s <- 4

        mbigm <- rbind(mbigm, c("M", tp, rownames(lodM.m)[x], rownames(lodM.m)[y], s))

        mcnt <- mcnt + 1

        l2 <- map[all[rownames(lodM.m)[y]], "cloc"]
       # if(type == 1) drawspline(locs[l1, 1:2], locs[l2, 1:2], col = rgb(1,0,0,0.5), lty=1, lwd=s)
       # if(type == 2) drawspline(locs[l1, 1:2], locs[l2, 1:2], col = rgb(0.5,0.5,0.5,0.5), lty=2, lwd=s)
       # if(type == 3) drawspline(locs[l1, 1:2], locs[l2, 1:2], col = rgb(0,0,1,0.5), lty=1, lwd=s)
        return
      }
    }
    #points(locs[l1, 1:2], pch = pch, col = col)
  }
  for(x in allV){
    l1 <- map[x, "cloc"]
    #points(locs[l1, 1:2], pch = 1, col = "black")
  }

  for(x in allS){
    l1 <- map[x, "cloc"]
    #points(locs[l1, 1:2], pch = 20, col = "green")
  }
  #legend("topleft", c("Vita:Vita", "Vita:Soma", "Soma:Soma"), lwd=2, lty=c(1,2,1), col = c(rgb(1,0,0,0.5), rgb(0.5,0.5,0.5,0.5), rgb(0,0,1,0.5)), bty="n")
  #legend("topright", c("Vita", "Soma"), pch=c(1,20), col = c("black", "green"), bty="n")

### FEMALES

  #plot(locs[,1:2], pch = 19, cex=0.1, col = locs[, "col"],xlab="", ylab="", xaxt="n", yaxt="n", bty="n", main = paste0("Female", tp))
  #text(1.05 * locs[chr.mids, c(1:2)], names(chr.mids))
  fcnt <- 0
  for(x in 1:(nrow(lodM.f)-1)){
    l1 <- map[all[rownames(lodM.f)[x]], "cloc"]
    pch <- 1; col = "black";
    if(grepl("Soma", rownames(lodM.f)[x])){ pch <- 20; col = "green";}

    for(y in (x+1):nrow(lodM.f)){
      m1 <- gsub("Soma", "", gsub("Vita", "", rownames(lodM.f)[x]))
      m2 <- gsub("Soma", "", gsub("Vita", "", rownames(lodM.f)[y]))
      m1 <- substr(m1, 1, nchar(m1)-1)
      m2 <- substr(m2, 1, nchar(m2)-1)

      type <- 1
      if(grepl("Soma", rownames(lodM.m)[x]) && grepl("Vita", rownames(lodM.m)[y])) type <- 2
      if(grepl("Vita", rownames(lodM.m)[x]) && grepl("Soma", rownames(lodM.m)[y])) type <- 2
      if(grepl("Soma", rownames(lodM.m)[x]) && grepl("Soma", rownames(lodM.m)[y])) type <- 3

      if(m1 != m2 && lodM.f[x,y] > 2.0){
        cat(all[rownames(lodM.f)[x]], " ",all[rownames(lodM.f)[y]], "\n")
        s <- 1
        fcnt <- fcnt + 1
        if(lodM.f[x,y] > 3.8 && lodM.f[x,y] <= 4.2) s <- 2
        if(lodM.f[x,y] > 4.2 && lodM.f[x,y] <= 4.5) s <- 3
        if(lodM.f[x,y] > 4.5 ) s <- 4

        mbigm <- rbind(mbigm, c("F", tp, rownames(lodM.f)[x], rownames(lodM.f)[y], s))

        l2 <- map[all[rownames(lodM.f)[y]], "cloc"]
        #if(type == 1) drawspline(locs[l1, 1:2], locs[l2, 1:2], col = rgb(1,0,0,0.5), lty=1, lwd=s)
        #if(type == 2) drawspline(locs[l1, 1:2], locs[l2, 1:2], col = rgb(0.5,0.5,0.5,0.5), lty=2, lwd=s)
        #if(type == 3) drawspline(locs[l1, 1:2], locs[l2, 1:2], col = rgb(0,0,1,0.5), lty=1, lwd=s)
        return
      }
    }
    #points(locs[l1, 1:2], pch = pch, col = col)
  }
  for(x in allV){
    l1 <- map[x, "cloc"]
    #points(locs[l1, 1:2], pch = 1, col = "black")
  }

  for(x in allS){
    l1 <- map[x, "cloc"]
    #points(locs[l1, 1:2], pch = 20, col = "green")
  }


  #legend("topright", c("Vita:Vita", "Vita:Soma", "Soma:Soma"), lwd=2, lty=c(1,2,1), col = c(rgb(1,0,0,0.5), rgb(0.5,0.5,0.5,0.5), rgb(0,0,1,0.5)), bty="n")
  #legend("topleft", c("Vita", "Soma"), pch=c(1,20), col = c("black", "green"), bty="n")

  #
  #legend("bottomleft", c("3.8 - 4.2", "4.2 - 4.5", "4.5+"), lwd=c(2,3,4), col = c("black"), bty="n")
  #legend("bottomright", c(paste0("Male=", mcnt), paste0("Female=", fcnt)), lwd=c(2,3,4), col = c("black"), bty="n")

#  dev.off()
}

mbigm <- data.frame(mbigm)

num <- c()
overlap <- c()
for(x in 1:4){
  num <- rbind(num , c("M",x, table(mbigm[mbigm[,5] >= x & mbigm[,1] == "M",2])[c("42", "365", "740", "905")]))
  num <- rbind(num , c("F",x, table(mbigm[mbigm[,5] >= x & mbigm[,1] == "F",2])[c("42", "365", "740", "905")]))

  for(tp in c("42", "365", "740", "905")){
    maxO <- min(as.numeric(num[which(num[,2] == as.character(x)), tp]))
    ml <- mbigm[mbigm[,5] >= x & mbigm[,1] == "M" & mbigm[,2] == tp,]
    mn <- paste0(ml[,3], ":", ml[,4])

    fl <- mbigm[mbigm[,5] >= x & mbigm[,1] == "F" & mbigm[,2] == tp,]
    fn <- paste0(fl[,3], ":", fl[,4])

    overlap <- rbind(overlap, c(tp, x, length(which(mn %in% fn)), maxO))
    cat("TP:", tp, "@", x, ", overlap=", length(which(mn %in% fn)),"/", maxO, ", olap=", paste0(mn[which(mn %in% fn)], sep = ", "),  "\n",sep="")
  }
}

write.table(num, "GxG_numSignificant.txt", sep = "\t", row.names=FALSE)
write.table(overlap, "GxG_overlapMvsF.txt", sep = "\t", row.names=FALSE)

for(x in 1:4){
  ml <- mbigm[mbigm[,5] >= x & mbigm[,1] == "M",]
  tp42m <- paste0(ml[which(ml[,2] == 42),3], ":", ml[which(ml[,2] == 42),4])
  tp365m <- paste0(ml[which(ml[,2] == 365),3], ":", ml[which(ml[,2] == 365),4])
  tp740m <- paste0(ml[which(ml[,2] == 740),3], ":", ml[which(ml[,2] == 740),4])
  tp905m <- paste0(ml[which(ml[,2] == 905),3], ":", ml[which(ml[,2] == 905),4])

  cat("M ", x, " 365 42 ", length(which(tp365m %in% tp42m)), " ", min(length(tp365m), length(tp42m)), "\n",sep="")
  cat("M ", x, " 740 42 ", length(which(tp740m %in% tp42m)), " ", min(length(tp740m), length(tp42m)), "\n",sep="")
  cat("M ", x, " 905 42 ", length(which(tp905m %in% tp42m)), " ", min(length(tp905m), length(tp42m)), "\n",sep="")
  cat("M ", x, " 740 365 ", length(which(tp740m %in% tp365m)), " ", min(length(tp740m), length(tp365m)), "\n",sep="")
  cat("M ", x, " 905 365 ", length(which(tp905m %in% tp365m)), " ", min(length(tp905m), length(tp365m)), "\n",sep="")
  cat("M ", x, " 905 740 ", length(which(tp905m %in% tp740m)), " ", min(length(tp905m), length(tp740m)), "\n",sep="")

  ml <- mbigm[mbigm[,5] >= x & mbigm[,1] == "F",]
  tp42m <- paste0(ml[which(ml[,2] == 42),3], ":", ml[which(ml[,2] == 42),4])
  tp365m <- paste0(ml[which(ml[,2] == 365),3], ":", ml[which(ml[,2] == 365),4])
  tp740m <- paste0(ml[which(ml[,2] == 740),3], ":", ml[which(ml[,2] == 740),4])
  tp905m <- paste0(ml[which(ml[,2] == 905),3], ":", ml[which(ml[,2] == 905),4])

  cat("F ", x, " 365 42 ", length(which(tp365m %in% tp42m)), " ", min(length(tp365m), length(tp42m)), "\n",sep="")
  cat("F ", x, " 740 42 ", length(which(tp740m %in% tp42m)), " ", min(length(tp740m), length(tp42m)), "\n",sep="")
  cat("F ", x, " 905 42 ", length(which(tp905m %in% tp42m)), " ", min(length(tp905m), length(tp42m)), "\n",sep="")
  cat("F ", x, " 740 365 ", length(which(tp740m %in% tp365m)), " ", min(length(tp740m), length(tp365m)), "\n",sep="")
  cat("F ", x, " 905 365 ", length(which(tp905m %in% tp365m)), " ", min(length(tp905m), length(tp365m)), "\n",sep="")
  cat("F ", x, " 905 740 ", length(which(tp905m %in% tp740m)), " ", min(length(tp905m), length(tp740m)), "\n",sep="")
}

