setwd("C:/Users/rqdt9/OneDrive - Northumbria University - Production Azure AD/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))


lods.cM <- read.table("progressiveMapping_all.txt", sep = "\t", check.names=FALSE)
lods.mM <- read.table("progressiveMapping_males.txt", sep = "\t", check.names=FALSE)
lods.fM <- read.table("progressiveMapping_females.txt", sep = "\t", check.names=FALSE)

lods.c <- apply(lods.cM,2,max)
lods.m <- apply(lods.mM,2,max)
lods.f <- apply(lods.fM,2,max)

# Plot the QTL profile
chrs <- unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",1))
positions <- as.numeric(unlist(lapply(strsplit(colnames(pull.geno(mcross)), "_"),"[",2)))
chr.l <- c()
for(chr in unique(chrs)){ chr.l <- c(chr.l, max(positions[which(chrs==chr)])); }
names(chr.l) <- unique(chrs)

plot(c(-1000000, 5000000 + sum(chr.l) + 19 * 25000000), c(0,12), main = "QTL profiles for longevity", t = 'n', xaxt='n', xlab="Chromosome", ylab="LOD", xaxs="i", yaxs="i", las=2)
#abline(h = seq(0,10,2), lty=2)
s <- 0
h <- c()
for(chr in unique(chrs)){
  rect(s, 0, s + chr.l[chr] + 5000000, 12, col = rgb(0.9,0.9,0.9,0.5),border=NA)
  points(positions[chrs==chr] + s, lods.c[chrs == chr], t = 'l', lwd=3, col="black")
  points(positions[chrs==chr] + s, lods.f[chrs == chr], t = 'l', lwd=2, col="pink")
  points(positions[chrs==chr] + s, lods.m[chrs == chr], t = 'l', lwd=2, col="blue")
  h <- c(h, s + (chr.l[chr]/2) + 25000000)
  s <- s + chr.l[chr] + 25000000
}
axis(1, at = h, unique(chrs))
abline(h = c(3.65, 4.25, 4.95), col = c("gold", "orange", "green"), lty=2, lwd=2)
legend("topright", c("All", "Only females", "Only males"), lwd=c(3,2,2), col=c("black", "pink", "blue"), bg="white")
legend("topleft", c("P(fdr) = 0.1", "P(bf) = 0.05", "P(bf) = 0.01"), lwd=1, lty=2, col=c("gold", "orange", "green"), bg="white")
box()

