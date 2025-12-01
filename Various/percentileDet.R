
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder")
regions <- read.csv("Vita_For_DA_28Nov2025.csv")
regions[,1] <- unlist(lapply(strsplit(regions[,1], " "), "[", 1))
regions[regions[,4] == "C",4] = "all"
regions[regions[,4] == "M",4] = "males"
regions[regions[,4] == "F",4] = "females"


setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
map <- read.table("genetic_map.txt", sep = "\t")


mtable <- c()
for(x in 1:nrow(regions)){
  ### Data input
  region.name <- regions[x,1]
  region.tp <- paste0("> ", regions[x,3])
  region.chr <- regions[x, 2]
  region.type <- regions[x, 4]


  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/regionsAug25/genes")
  genes <- read.table(paste0("protein_coding_genes_",region.name,".txt"), sep = "\t",  quote = "")
  pos <- genes[, "start_position"] + (genes[, "end_position"] - genes[, "start_position"])

  setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
  lods <- read.table(paste0("progressiveMapping_",region.type,"20D.txt"), sep = "\t", check.names=FALSE)

  onChr <- rownames(map)[which(map[, 1] == region.chr)]
  peak <- max(lods[region.tp, onChr])

  for(p in 1:length(pos)){
    n <-  genes[p, "external_gene_name"]
    m1 <- onChr[max(which(map[onChr,2] < pos[p]))]
    m2 <- onChr[min(which(map[onChr,2] > pos[p]))]
    lodB <- lods[region.tp, m1]
    lodA <- lods[region.tp, m2]
    lodDrop <- peak - mean(c(lodA,lodB),na.omit=TRUE)
    cat(n, " ",  pos[p], " ", m1, " ", m2, " ", lodB, " ", lodA, " ", lodDrop, "\n")
    mtable <- rbind(mtable, c(region = region.name, name = n, position = pos[p], m1 = m1, m2 = m2, lodDrop = round(lodDrop,2)))
  }
}

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder")
write.table(mtable, "lodDrops.txt", sep = "\t", quote = FALSE, row.names = FALSE)


############ SOma


setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder")
regions <- read.csv("Soma_For_DA_28Nov2025.csv", sep = "\t")
regions[,1] <- unlist(lapply(strsplit(regions[,1], " "), "[", 1))

setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
map <- read.table("genetic_map.txt", sep = "\t")


mtable <- c()
for(x in 1:nrow(regions)){
  ### Data input
  region.name <- regions[x,1]
  region.tp <- regions[x,3]
  region.chr <- regions[x, 2]
  region.type <- regions[x, 4]
  region.bw <- regions[x, 5]


  setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/regionsAug25/genes")
  genes <- read.table(paste0("protein_coding_genes_",region.name,".txt"), sep = "\t",  quote = "")
  pos <- genes[, "start_position"] + (genes[, "end_position"] - genes[, "start_position"])

  setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")
  lods <- t(read.table(paste0("CTL_BW",region.bw,"_T",region.tp,"_",region.type,".txt"), sep = "\t", check.names=FALSE))


  onChr <- rownames(map)[which(map[, 1] == region.chr)]
  onChr <-  onChr[which(onChr %in% colnames(lods))]
  peak <- max(lods[1, onChr])

  for(p in 1:length(pos)){
    n <-  genes[p, "external_gene_name"]
    m1 <- onChr[max(which(map[onChr,2] < pos[p]))]
    m2 <- onChr[min(which(map[onChr,2] > pos[p]))]
    lodB <- NA
    lodA <- NA
    if(!is.na(m1)) lodB <- lods[1, m1]
    if(!is.na(m2)) lodA <- lods[1, m2]
    lodDrop <- peak - mean(c(lodA,lodB),na.rm=TRUE)
    cat(n, " ",  pos[p], " ", m1, " ", m2, " ", lodB, " ", lodA, " ", lodDrop, "\n")
    mtable <- rbind(mtable, c(region = region.name, name = n, position = pos[p], m1 = m1, m2 = m2, lodDrop = round(lodDrop,2)))
  }
}

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder")
write.table(mtable, "lodDrops_SOma.txt", sep = "\t", quote = FALSE, row.names = FALSE)

