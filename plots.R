setwd("C:/Users/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/merged")

gts4way <- read.table("gts4way.txt", sep="\t")
gts4wayRqtl <- read.table("gts4way.rqtl.txt", sep="\t")
map <- read.table("map.gts4way.txt", sep="\t")
ind <- read.table("ind.gts4way.txt", sep="\t")

plot(c(1,19), c(0, max(map[, "Position"])), t = 'n', xlab="Chromosome", ylab="Position (Mbp)",xaxt='n', yaxt='n')
for(x in 1:19){
  points(c(x,x), c(0, max(map[which(map[,2] == x),3])), t = 'l', lwd=2)
}
apply(map, 1, function(x){
  iX <- which(c(1:19) == x[2])
  if(length(iX) > 0){
    off <- ((as.numeric(grepl("Pat", x[6])) - 0.5) / 5)
    points(iX+off, x[3], pch="-", col=c("blue", "pink")[as.numeric(grepl("Mat", x[6])) + 1], cex=2)
  }
})
axis(1, at = 1:19, 1:19, las=1)
axis(2, at = seq(0, max(map[, "Position"]), 25000000), seq(0, max(map[, "Position"]), 25000000) / 1000000, las=1)


