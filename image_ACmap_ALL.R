library(RColorBrewer)
library(svglite)

setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
map <- read.table("genetic_map.txt", sep = "\t")

lods.m.All <- read.table("progressiveMapping_males.txt", sep = "\t", check.names=FALSE)
#lods.m.All <- lods.m.All[c(TRUE,FALSE,FALSE), ]
lods.f.All <- read.table("progressiveMapping_females.txt", sep = "\t", check.names=FALSE)
#lods.f.All <- lods.f.All[c(TRUE,FALSE,FALSE), ]
lods.c.All <- read.table("progressiveMapping_all.txt", sep = "\t", check.names=FALSE)
#lods.c.All <- lods.c.All[c(TRUE,FALSE,FALSE), ]

threshold <- 2.5
# Do some basic checks
head(map)
rownames(map)[1:10]
colnames(map)

# Look at the chromosome column and count the number of markers on each
table(map[, "Chr"])

y.max <- max(map[, "Pos"])
y.min <- 0
x.min <- 1
x.max <- length(table(map[, "Chr"]))

plot(x = c(x.min, x.max),
     y = c(y.min, y.max),
     type = "n",
     xlab = "Chromosomes",
     ylab = "Position",
     main = "Genetic map UM-HET3",
     xaxt = "n",
     yaxt = "n"
     )
# Plot our own X axis, use the unique entries in the map matrix
axis(1, at = 1:x.max, unique(map[, "Chr"]))

y.labels <- seq(0, y.max, 10000000) / 1000000
axis(2, at = seq(0, y.max, 10000000), y.labels, las = 2)

for(i in 1:nrow(map)) {
  pos.x <- map[i, "Chr"]
  if(pos.x == "X") pos.x <- x.max
  points(x = pos.x, y = map[i, "Pos"], pch = "-")
}

chrs <- c(1:19, "X")
gap <- 10000000
chrs.length <- unlist(lapply(chrs, function(x){ max(map[map[, "Chr"] == x,"Pos"]) }))
names(chrs.length) <- chrs
l.x <- sum(chrs.length) + gap * (length(chrs)-1)


colz.c <- c("white", brewer.pal(9, "Greens")[-c(1:4)])
colz.m <- c("white", brewer.pal(9, "Blues")[-c(1:4)])
colz.f <- c("white", brewer.pal(9, "PuRd")[-c(1:4)])

thresholds <- c(0, 2, 3.65, 4.25, 4.95, 5.95, 100)
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
pdf(paste0("Figure_1_actuaryQTL_all.pdf"), width = 24, height = 12)


plot(c(1, l.x), c(1, 1+nrow(lods.c.All)), t = 'n', xlab = "", xaxt='n', ylab="", yaxt='n', xaxs="i", yaxs="i")
abline(h = 1:nrow(lods.c.All), col = rgb(0.5,0.5,0.5))
chr.s <- 0
chr.pe <- 0
c.i <- 1
for(chr in chrs){
  rect(chr.s, 0, chr.s + chrs.length[chr], 1+nrow(lods.c.All), col = c(rgb(1,1,1,0.5), rgb(1,1,1,0.5))[1 + (c.i %% 2)], border=NA); c.i <- c.i + 1;
  axis(1, at = chr.s + chrs.length[chr] / 2, chr)
  rect(chr.pe, 0, chr.s, 1+nrow(lods.c.All), col = "white", border=NA)
  abline(v = chr.s)
  abline(v = chr.s + chrs.length[chr])
  rr <- rownames(map)[which(map[, "Chr"] == chr)]
  rr.pos <- map[rr, "Pos"] + chr.s
  
  for(tp in 1:nrow(lods.c.All)){
    ii <- which(lods.c.All[tp,rr] > threshold)
    for(i in ii){
      if(i == 1){ x.s <- chr.s; }else{ x.s <- (rr.pos[i] + rr.pos[i-1])/2; }
      if(i == length(rr)){ x.e <- chr.s + chrs.length[chr]; }else{ x.e <- (rr.pos[i] + rr.pos[i+1])/2; }
      lod <- as.numeric(lods.c.All[tp,rr][i])
      #cat(x.s," ",  x.e, " ", lod, "\n")
      cN <- max(which(thresholds < lod))
      rect(x.s, tp+0.4, x.e, tp+0.6, col = colz.c[cN], border = NA)
    }
  }

  for(tp in 1:nrow(lods.f.All)){
    ii <- which(lods.f.All[tp,rr] > threshold)
    for(i in ii){
      if(i == 1){ x.s <- chr.s; }else{ x.s <- (rr.pos[i] + rr.pos[i-1])/2; }
      if(i == length(rr)){ x.e <- chr.s + chrs.length[chr]; }else{ x.e <- (rr.pos[i] + rr.pos[i+1])/2; }
      lod <- as.numeric(lods.f.All[tp,rr][i])
      #cat(x.s," ",  x.e, " ", lod, "\n")
      cN <- max(which(thresholds < lod))
      rect(x.s, tp+0.1, x.e, tp+0.3, col = colz.f[cN], border = NA)
    }
  }

  for(tp in 1:nrow(lods.m.All)){
    ii <- which(lods.m.All[tp,rr] > threshold)
    for(i in ii){
      if(i == 1){ x.s <- chr.s; }else{ x.s <- (rr.pos[i] + rr.pos[i-1])/2; }
      if(i == length(rr)){ x.e <- chr.s + chrs.length[chr]; }else{ x.e <- (rr.pos[i] + rr.pos[i+1])/2; }
      lod <- as.numeric(lods.m.All[tp,rr][i])
      #cat(x.s," ",  x.e, " ", lod, "\n")
      cN <- max(which(thresholds < lod))
      rect(x.s, tp+0.7, x.e, tp+0.90, col = colz.m[cN], border = NA)
    }
  }
  chr.pe <- chr.s + chrs.length[chr]
  chr.s <- chr.s + chrs.length[chr] + gap
}
axis(2, at = 0.5 + 1:nrow(lods.c.All), rownames(lods.c.All), las=2)
legend("topright", fill = c(colz.m[3],colz.c[3],colz.f[3]), c("Males", "Combined", "Females"), bg = "white")

dev.off()


