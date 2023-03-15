setwd("C:/Users/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3")

regions <- read.table("regions_4way_merged.txt", sep="\t", header=FALSE, row.names=1)
colnames(regions) <- c("Chr", "Proximal", "Distal")

library(biomaRt)

mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org")

mlist <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){
  r <- paste0(regions[x, "Chr"], ":",regions[x, "Proximal"],":", regions[x, "Distal"])
  mlist[[x]] <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", 
                                     "chromosome_name", "start_position", "end_position", "strand", "description"), 
                      filters = "chromosomal_region", values = r, mart = mart)
}
names(mlist) <- rownames(regions)

types <- unique(unlist(lapply(mlist, function(x){names(table(x[, "gene_biotype"]))})))
mm <- matrix(0, nrow(regions), length(types), dimnames=list(rownames(regions), types))

for(x in 1:nrow(regions)){
  tbl <- table(mlist[[x]][, "gene_biotype"])
  mm[x, names(tbl)] <- tbl
}
write.table(mm, "FeaturesVitaLoci.txt",sep="\t", quote=FALSE)

for(x in 1:nrow(regions)){
  mVEP <- read.csv(paste0("SNPVEPFeb23/",rownames(regions)[x], ".txt.vep"),sep="\t", skip=80, header=FALSE)
  
  mGenAge <- read.csv(paste0("genes/all/",rownames(regions)[x], ".txt"), sep="\t", header=TRUE, row.names = 1)
  
  colnames(mVEP) <- c("Uploaded_variation","Location","Allele","Gene","Feature","Feature_type","Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","Extra")
  isInteresting <- c(grep("HIGH", mVEP[,"Extra"]), grep("MODERATE", mVEP[,"Extra"]))
  genes <- unique(mVEP[isInteresting, "Gene"])
  mSub <- mVEP[isInteresting,]

  mm <- matrix(NA, length(genes), 13, dimnames = list(genes, c("external_gene_name", "isHIGH", "inGenAge", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "description", "genAge", "nLocations", "nFeatures", "variations")))
  for(gene in genes){
    inmList <- which(mlist[[x]][,1] == gene)
    if(length(inmList) == 1){
      mm[gene, "external_gene_name"] <- mlist[[x]][inmList,"external_gene_name"]
      mm[gene, "chromosome_name"] <- mlist[[x]][inmList,"chromosome_name"]
      mm[gene, "start_position"] <- mlist[[x]][inmList,"start_position"]
      mm[gene, "end_position"] <- mlist[[x]][inmList,"end_position"]
      mm[gene, "strand"] <- mlist[[x]][inmList,"strand"]
      mm[gene, "gene_biotype"] <- mlist[[x]][inmList,"gene_biotype"]
      mm[gene, "description"] <- mlist[[x]][inmList,"description"]
    }
    mm[gene, "nLocations"] <- length(unique(mSub[which(mSub[, "Gene"] == gene),"Location"]))
    mm[gene, "nFeatures"] <- length(unique(mSub[which(mSub[, "Gene"] == gene),"Feature"]))
    mm[gene, "variations"] <- paste0(mSub[which(mSub[, "Gene"] == gene),"Consequence"], collapse=";")
    
    isHIGH <- any(grepl("HIGH", mSub[which(mSub[, "Gene"] == gene),"Extra"]))
    if(isHIGH){
      mm[gene, "isHIGH"] <- "YES"
    }else{
      mm[gene, "isHIGH"] <- "NO"
    }
    
    if(gene %in% rownames(mGenAge)){
      mm[gene, "inGenAge"] <- "YES"
      mm[gene, "genAge"] <- mGenAge[gene, "genAge"]
    }else{
      mm[gene, "inGenAge"] <- "NO"
    }
  }
  
  # Filter some (Riken, Gm, is.NA)
  hasRIKEN <- grep("RIKEN cDNA", mm[, "description"])
  if(length(hasRIKEN) > 0){ mm <- mm[-hasRIKEN,] }
  isPredicted <- grep("predicted gene", mm[, "description"])
  if(length(isPredicted) > 0){ mm <- mm[-isPredicted,] }
  isNA <- which(is.na(mm[, "description"]))
  if(length(isNA) > 0){ mm <- mm[-isNA,] }
  
  colnames(mm)[1] <- "name"
  colnames(mm)[4] <- "chr"
  colnames(mm)[5] <- "start"
  colnames(mm)[6] <- "stop"
  write.table(mm, file = paste0("SNPVEPFeb23/",rownames(regions)[x], ".summary.txt"), sep = "\t", quote = FALSE,na="")
}
