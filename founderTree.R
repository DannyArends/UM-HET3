setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files/merged")
fvcf <- read.table("fvcfAll.txt", sep="\t",colClasses="character", header=TRUE)

gts <- apply(fvcf[,6:9],1,function(x){as.numeric(factor(x, levels = c("0/0", "0/1", "1/1")))})-2
rownames(gts) <- colnames(fvcf)[6:9]
plot(hclust(dist(gts)))
