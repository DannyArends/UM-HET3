#
# circle.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Code to generate a circle plot with splines showing GxG interactions between (and among) for Vita and Soma loci
#

library(qtl)

source("ActuarialMapping/adjustXprobs.R")

# Read cross object
mcross <- read.cross(format="csvr", file="DataSet/um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)
mcross <- adjustXprobs(mcross)

gtsp <- pull.genoprob(mcross)

lods.all <- read.table("DataSet/output/progressiveMapping_all20D.txt", sep = "\t", check.names = FALSE)
lods.f <- read.table("DataSet/output/progressiveMapping_females20D.txt", sep = "\t", check.names = FALSE)
lods.m <- read.table("DataSet/output/progressiveMapping_males20D.txt", sep = "\t", check.names = FALSE)

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

tp <- 42
lodM.m <- read.table(paste0("DataSet/output/vita_soma_interactions_2way_males_tp",tp,".txt"), sep = "\t")
for(x in 1:ncol(lodM.m)){
for(y in 1:x){
    lodM.m[y,x] <- lodM.m[x,y]
  }
} 
lodM.m <- lodM.m[names(all), names(all)]

lodM.f <- read.table(paste0("DataSet/output/vita_soma_interactions_2way_females_tp",tp,".txt"), sep = "\t")
for(x in 1:ncol(lodM.f)){
for(y in 1:x){
    lodM.f[y,x] <- lodM.f[x,y]
  }
} 
lodM.f <- lodM.f[names(all), names(all)]

map <- read.table("DataSet/output/genetic_map.txt", sep = "\t")
map <- cbind(map, cloc = NA)
map <- map[-c(893,892),]

toPos <- function(chr, pos){ return(chr.starts[chr] + round(pos / 1e6)) }

for(x in 1:nrow(map)){
  map[x, "cloc"] = toPos(map[x, "Chr"], map[x, "Pos"])
}

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

plot(locs[,1:2], pch = 19, cex=0.1, col = locs[, "col"],xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
maxS <- max(lods.all[1,])
off <- 0.99
for(y in seq(1, nrow(lods.all), 1)){
  for(x in 1:nrow(map)){
    p <- locs[map[x, "cloc"], 1:2]
    lod <- lods.all[y,x]
    if(lod > 3){
      points(off * p, col = "green", pch = 19, cex=0.2)
    }
  }
  off <- off - 0.0025
}


off <- off - 0.01
foff <- off
points(off * locs[,1:2], pch = 19, cex=0.1, col = locs[, "col"])
off <- off - 0.01
for(y in seq(1, nrow(lods.f), 1)){
  for(x in 1:nrow(map)){
    p <- locs[map[x, "cloc"], 1:2]
    lod <- lods.f[y,x]
    if(lod > 3){
      points(off * p, col = "red", pch = 19, cex=0.2)
    }
  }
  off <- off - 0.0025
}
off <- off - 0.01

moff <- off
points(off * locs[,1:2], pch = 19, cex=0.1, col = locs[, "col"])
off <- off - 0.01
for(y in seq(1, nrow(lods.m), 1)){
  for(x in 1:nrow(map)){
    p <- locs[map[x, "cloc"], 1:2]
    lod <- lods.m[y,x]
    if(lod > 3){
      points(off * p, col = "blue", pch = 19, cex=0.2)
    }
  }
  off <- off - 0.0025
}
off <- off - 0.01
points(off * locs[,1:2], pch = 19, cex=0.1, col = locs[, "col"])

for(x in 1:(nrow(lodM.m)-1)){
  for(y in (x+1):nrow(lodM.m)){
    m1 <- gsub("Soma", "", gsub("Vita", "", rownames(lodM.m)[x]))
    m2 <- gsub("Soma", "", gsub("Vita", "", rownames(lodM.m)[y]))
    m1 <- substr(m1, 1, nchar(m1)-1)
    m2 <- substr(m2, 1, nchar(m2)-1)

    if(m1 != m2 && lodM.m[x,y] > 5){
      cat(all[rownames(lodM.m)[x]], " ",all[rownames(lodM.m)[y]], "\n")
      l1 <- map[all[rownames(lodM.m)[x]], "cloc"]
      l2 <- map[all[rownames(lodM.m)[y]], "cloc"]
      drawspline(moff *locs[l1, 1:2],moff * locs[l2, 1:2], col = rgb(0,0,1,0.5))
      return
    }
  }
}

for(x in 1:(nrow(lodM.f)-1)){
  for(y in (x+1):nrow(lodM.f)){
    m1 <- gsub("Soma", "", gsub("Vita", "", rownames(lodM.f)[x]))
    m2 <- gsub("Soma", "", gsub("Vita", "", rownames(lodM.f)[y]))
    m1 <- substr(m1, 1, nchar(m1)-1)
    m2 <- substr(m2, 1, nchar(m2)-1)

    if(m1 != m2 && lodM.f[x,y] > 5){
      cat(all[rownames(lodM.f)[x]], " ",all[rownames(lodM.f)[y]], "\n")
      l1 <- map[all[rownames(lodM.f)[x]], "cloc"]
      l2 <- map[all[rownames(lodM.f)[y]], "cloc"]
      drawspline(foff * locs[l1, 1:2], foff * locs[l2, 1:2], col = rgb(1,0,0,0.5))
      return
    }
  }
}



