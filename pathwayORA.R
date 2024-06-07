setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3")

regions <- read.table("regions_4way_merged_May24.txt", sep="\t", header=FALSE, row.names=1)
colnames(regions) <- c("Chr", "Proximal", "Distal", "Old")

protein_coding <- vector("list", nrow(regions))
prioritized <- vector("list", nrow(regions))

bg <- c()
prio <- c()
prioH <- c()
for(x in 1:nrow(regions)){
  protein_coding[[x]] <- read.csv(paste0("May2024/genes/protein_coding_genes_", rownames(regions)[x], ".txt"), sep = "\t")
  prioritized[[x]] <- read.csv(paste0("May2024/summary/", rownames(regions)[x], ".summary.txt"), sep = "\t")

  iiHigh <- which(prioritized[[x]][,3] == "HIGH")
  iiMod <- which(prioritized[[x]][,3] != "NONE")

  bg <- c(bg, protein_coding[[x]][,1])
  prio <- c(prio, prioritized[[x]][iiMod,1])
  prioH <- c(prioH, prioritized[[x]][iiHigh,1])
}

cat(unique(bg), file = "ORAbackground.txt",sep="\n", collapse = "\n")
cat(unique(prio), file = "ORAprio.txt",sep="\n", collapse = "\n")
cat(unique(prioH), file = "ORAprioH.txt",sep="\n", collapse = "\n")

regions <- read.table("regions_4way_merged_March23.txt", sep="\t", header=FALSE, row.names=1)
colnames(regions) <- c("Chr", "Proximal", "Distal", "Old")
regions <- regions[which(regions[,4] == 1),]

bg <- c()
prio <- c()
prioH <- c()
for(x in 1:nrow(regions)){
  protein_coding[[x]] <- read.csv(paste0("May2024/genes/protein_coding_genes_", rownames(regions)[x], ".txt"), sep = "\t")
  prioritized[[x]] <- read.csv(paste0("May2024/summary/", rownames(regions)[x], ".summary.txt"), sep = "\t")

  iiHigh <- which(prioritized[[x]][,3] == "HIGH")
  iiMod <- which(prioritized[[x]][,3] != "NONE")

  bg <- c(bg, protein_coding[[x]][,1])
  prio <- c(prio, prioritized[[x]][iiMod,1])
  prioH <- c(prioH, prioritized[[x]][iiHigh,1])
}

cat(unique(bg), file = "ORAbackground_OLD.txt",sep="\n", collapse = "\n")
cat(unique(prio), file = "ORAprio_OLD.txt",sep="\n", collapse = "\n")
cat(unique(prioH), file = "ORAprioH_OLD.txt",sep="\n", collapse = "\n")


regions <- read.table("regions_4way_merged_March23.txt", sep="\t", header=FALSE, row.names=1)
colnames(regions) <- c("Chr", "Proximal", "Distal", "Old")
regions <- regions[which(regions[,4] == 0),]

bg <- c()
prio <- c()
prioH <- c()
for(x in 1:nrow(regions)){
  protein_coding[[x]] <- read.csv(paste0("genes/protein_coding_genes_", rownames(regions)[x], ".txt"), sep = "\t")
  prioritized[[x]] <- read.csv(paste0("RegionsMarch14/", rownames(regions)[x], ".summary.txt"), sep = "\t")

  iix <- which(prioritized[[x]][,3] == "HIGH")

  bg <- c(bg, protein_coding[[x]][,1])
  prio <- c(prio, prioritized[[x]][,1])
  prioH <- c(prioH, prioritized[[x]][iix,1])
}

cat(unique(bg), file = "ORAbackground_Young.txt",sep="\n", collapse = "\n")
cat(unique(prio), file = "ORAprio_Young.txt",sep="\n", collapse = "\n")
cat(unique(prioH), file = "ORAprioH_Young.txt",sep="\n", collapse = "\n")

