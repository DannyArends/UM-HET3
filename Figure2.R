library(RColorBrewer)
library(svglite)

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
lods.cI <- read.table("progressiveMapping_INT.txt", sep = "\t",check.names=FALSE)
lods.cI <- lods.cI[c(TRUE,FALSE,FALSE), ]

map <- read.table("genetic_map.txt", sep = "\t")
threshold <- 2.5

chrs <- c(1:19, "X")
gap <- 10000000
chrs.length <- as.numeric(unlist(lapply(chrs, function(x){ max(map[map[, "Chr"] == x,"Pos"]) })))
names(chrs.length) <- chrs
l.x <- sum(chrs.length) + gap * (length(chrs)-1)

y.max <- max(map[, "Pos"])
y.min <- 0
x.min <- 1
x.max <- length(table(map[, "Chr"]))

colz.c <- c("white", brewer.pal(9, "Greens")[-c(1:4)])

thresholds <- c(0, 2, 3.65, 4.25, 4.95, 5.95, 100)

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/0000_ITP_BioRxiv_Tables_Files/Figures")
svglite(paste0("Figure_2_Haplo_x_Sex.svg"), width = 24, height = 12)

plot(c(1, l.x), c(1, 1+nrow(lods.cI)), t = 'n', xlab = "", xaxt='n', ylab="", yaxt='n', xaxs="i", yaxs="i", main = "Haplotype x Sex interaction")
abline(h = 1:nrow(lods.cI))
chr.s <- 0
chr.pe <- 0
c.i <- 1
for(chr in chrs){
  rect(chr.s, 0, chr.s + chrs.length[chr], 1+nrow(lods.cI), col = c(rgb(1,1,1,0.5), rgb(0.9,0.9,0.9,0.5))[1 + (c.i %% 2)]); c.i <- c.i + 1;
  axis(1, at = chr.s + chrs.length[chr] / 2, chr)
  rect(chr.pe, 0, chr.s, 1+nrow(lods.cI), col = "white")
  abline(v = chr.s)
  abline(v = chr.s + chrs.length[chr])
  rr <- rownames(map)[which(map[, "Chr"] == chr)]
  rr.pos <- as.numeric(map[rr, "Pos"]) + chr.s
  
  for(tp in 1:nrow(lods.cI)){
    ii <- which(lods.cI[tp,rr] > threshold)
    for(i in ii){
      if(i == 1){ x.s <- chr.s; }else{ x.s <- (rr.pos[i] + rr.pos[i-1])/2; }
      if(i == length(rr)){ x.e <- chr.s + chrs.length[chr]; }else{ x.e <- (rr.pos[i] + rr.pos[i+1])/2; }
      lod <- as.numeric(lods.cI[tp,rr][i])
      #cat(x.s," ",  x.e, " ", lod, "\n")
      cN <- max(which(thresholds < lod))
      rect(x.s, tp, x.e, tp+1.0, col = colz.c[cN], border = NA)
    }
  }

  chr.pe <- chr.s + chrs.length[chr]
  chr.s <- chr.s + chrs.length[chr] + gap
}
axis(2, at = 0.5 + 1:nrow(lods.cI), rownames(lods.cI), las=2)
dev.off()



