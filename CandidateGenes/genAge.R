#
# genAge.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Mapping our mouse genes to see if any GenAge homologous genes in C. elegans, Human, D. melanogaster, and S. cerevisiae are known
#

library(biomaRt)

GAmodel <- read.table("DataSet/genage/genage_models.csv",sep=",", header=TRUE)
GAhuman <- read.table("DataSet/genage/genage_human.csv",sep=",", header=TRUE)

regions <- read.table("DataSet/regions/regions_4way_merged.txt", sep="\t", header=FALSE, row.names=1)
colnames(regions) <- c("Chr", "Proximal", "Distal")

mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org")
for(x in 1:nrow(regions)){
  r <- paste0(regions[x, "Chr"],":", regions[x, "Proximal"], ":", regions[x, "Distal"])
  ortho <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_associated_gene_name", "celegans_homolog_associated_gene_name", "dmelanogaster_homolog_associated_gene_name", "scerevisiae_homolog_associated_gene_name"), filters = "chromosomal_region", values = r, mart = mart)
  descr <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "gene_biotype", "mgi_id", "mgi_symbol", "description", "mgi_description"), filters = "chromosomal_region", values = r, mart = mart)

  isGM <- grep("^Gm", descr[,"external_gene_name"])
  if(length(isGM) > 0) descr <- descr[-isGM,]

  isRIKEN <- grep("^RIKEN", descr[,"mgi_description"])
  if(length(isRIKEN) > 0) descr <- descr[-isRIKEN,]

  isNA <- which(is.na(descr[, "description"]))
  if(length(isNA) > 0){ descr <- descr[-isNA,] }

  write.table(descr, paste0("genes/all/",rownames(regions)[x], "_all.txt"), sep = "\t", quote = FALSE, na = "", row.names=FALSE)

  inGenAge <- c()
  for(co in c(4,5,6)){
    species <- gsub("_homolog_associated_gene_name","", colnames(ortho)[co])
    mHasG <- unique(ortho[which(ortho[, co] %in% GAmodel[,2]),1])
    if(length(mHasG) > 0){
      inGenAge <- rbind(inGenAge, cbind(species, mHasG))
    }
  }
  hHasG <- unique(ortho[which(ortho[, 3] %in% GAhuman[,2]),1])
  if(length(hHasG) > 0){
    inGenAge <- rbind(inGenAge, cbind("Human", hHasG))
  }

  descr <- cbind(descr, genAge = NA)
  for(ge in unique(inGenAge[,2])){
    specV <- paste0(as.character(inGenAge[which(inGenAge[,2] == ge),1]), collapse = ";")
    descr[which(descr[,"ensembl_gene_id"] == ge), "genAge"] <- specV
  }
  
  ## Make sure the "genes/all" folder exists ("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/2025")
  rname <- paste0("DataSet/output/genes/all/",rownames(regions)[x], "_genAge.txt")
  nFeat <- nrow(descr)
  nAll <- sum(descr[, "gene_biotype"] == "protein_coding")
  
  descr <- descr[which(!is.na(descr[, "genAge"])),]
  write.table(descr, rname, sep = "\t", quote = FALSE, na = "", row.names=FALSE)
  cat("Done",rownames(regions)[x], "in GenAge:", length(unique(inGenAge[,2])), "/", nAll, "total Features:", nFeat, "\n")
}

#Done Vita1a in GenAge: 1 / 14 total Features: 16 
#Done Vita1b in GenAge: 8 / 80 total Features: 97 
#Done Vita1c in GenAge: 14 / 289 total Features: 330 
#Done Vita2a in GenAge: 16 / 377 total Features: 479 
#Done Vita2b in GenAge: 18 / 246 total Features: 291 
#Done Vita2c in GenAge: 42 / 544 total Features: 632 
#Done Vita3a in GenAge: 43 / 628 total Features: 714 
#Done Vita3b in GenAge: 38 / 532 total Features: 596 
#Done Vita4a in GenAge: 16 / 262 total Features: 331 
#Done Vita5a in GenAge: 14 / 213 total Features: 252 
#Done Vita6a in GenAge: 13 / 143 total Features: 162 
#Done Vita6b in GenAge: 12 / 121 total Features: 144 
#Done Vita9a in GenAge: 14 / 278 total Features: 328 
#Done Vita9b in GenAge: 15 / 259 total Features: 291 
#Done Vita9c in GenAge: 10 / 206 total Features: 234 
#Done Vita10a in GenAge: 48 / 421 total Features: 475 
#Done Vita11a in GenAge: 17 / 154 total Features: 186 
#Done Vita11b in GenAge: 92 / 896 total Features: 1071 
#Done Vita12a in GenAge: 12 / 172 total Features: 497 
#Done Vita13a in GenAge: 14 / 119 total Features: 135 
#Done Vita14a in GenAge: 4 / 61 total Features: 86 
#Done Vita15a in GenAge: 28 / 307 total Features: 367 
#Done Vita17a in GenAge: 58 / 721 total Features: 878 
#Done Vita18a in GenAge: 17 / 178 total Features: 213 
#Done VitaXa in GenAge: 22 / 313 total Features: 422 
#Done VitaXb in GenAge: 11 / 111 total Features: 132 



