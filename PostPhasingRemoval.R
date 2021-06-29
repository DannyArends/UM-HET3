setwd("C:/Github/UM-HET3/files")
gts4way.fill <- read.table("merged/gts4way.filled.GN2.txt", sep="\t")
gts4way.rqtl <- read.table("merged/gts4way.rqtl.mai2021.txt", sep="\t")

# Use correlation to catch weird markers
mcor <- cor(t(gts4way.fill[,-c(1:2)]), use="pair")
badM <- names(which(apply(apply(mcor,1,is.na),1,sum) > 50)) # Some NAs are fine, but more than 50 means you're a bad marker

gts4way.rqtl <- gts4way.rqtl[-which(rownames(gts4way.rqtl) %in% badM),]
gts4way.fill <- gts4way.fill[-which(rownames(gts4way.fill) %in% badM),]

mcor <- cor(t(gts4way.fill[,-c(1:2)]), use="pair")

chr <- 6
image(which(gts4way.fill[,1] == chr), which(gts4way.fill[,1] == chr), mcor[gts4way.fill[,1] == chr, gts4way.fill[,1] == chr])
box()

# Chromsome 6 weird 434:437
# Chromsome 7 weird 494:497
# Chromsome 12 weird 849:850
# Chromsome 17 weird 1061:1065

ii <- c(434:437, 494:497, 847:850, 1061:1065)
gts4way.rqtl <- gts4way.rqtl[-ii,]
gts4way.fill <- gts4way.fill[-ii,]

mcor <- cor(t(gts4way.fill[,-c(1:2)]), use="pair")

image(1:nrow(mcor), 1:ncol(mcor), mcor)


chr <- 12
image(which(gts4way.fill[,1] == chr), which(gts4way.fill[,1] == chr), mcor[gts4way.fill[,1] == chr, gts4way.fill[,1] == chr])
box()

# too much missing data
tooMM <- names(which(apply(apply(gts4way.fill[-c(1,2),],1,is.na),2,sum)/ncol(gts4way.fill) > 0.7))

gts4way.rqtl <- gts4way.rqtl[-which(rownames(gts4way.rqtl) %in% tooMM),]
gts4way.fill <- gts4way.fill[-which(rownames(gts4way.fill) %in% tooMM),]

mcor <- cor(t(gts4way.fill[,-c(1:2)]), use="pair")
image(1:nrow(mcor), 1:ncol(mcor), mcor)

write.table(gts4way.rqtl, "merged/gts4way.filled.cleaned.GN2.txt", quote=FALSE, sep="\t", na = "")

map <- read.table("merged/map.gts4way.mai2021.txt", sep="\t")
map <- map[rownames(gts4way.rqtl),]
write.table(map, "merged/map.gts4way.filled.cleaned.GN2.txt", quote=FALSE, sep="\t", na = "")
