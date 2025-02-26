setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3")

regions <- read.table("regions_4way_2025.txt", sep="\t", header=FALSE, row.names=1)
colnames(regions) <- c("Chr", "Proximal", "Top", "Distal")


library(biomaRt)

mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org")

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/2025")

mlist <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){
  r <- paste0(regions[x, "Chr"], ":",regions[x, "Proximal"],":", regions[x, "Distal"])
  mlist[[x]] <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", 
                                     "chromosome_name", "start_position", "end_position", "strand", "description", "mgi_description"), 
                      filters = "chromosomal_region", values = r, mart = mart)
  isGM <- grep("^Gm", mlist[[x]][,"external_gene_name"])
  if(length(isGM) > 0) mlist[[x]] <- mlist[[x]][-isGM,]

  isRIKEN <- grep("^RIKEN", mlist[[x]][,"mgi_description"])
  if(length(isRIKEN) > 0) mlist[[x]] <- mlist[[x]][-isRIKEN,]

  isNA <- which(is.na(mlist[[x]][, "description"]))
  if(length(isNA) > 0){ mlist[[x]] <- mlist[[x]][-isNA,] }

  write.table(mlist[[x]][which(mlist[[x]][, "gene_biotype"] == "protein_coding"),], 
              paste0("genes/protein_coding_genes_", rownames(regions)[x], ".txt"),sep="\t",quote=FALSE)
}
names(mlist) <- rownames(regions)

types <- unique(unlist(lapply(mlist, function(x){names(table(x[, "gene_biotype"]))})))
mm <- matrix(0, nrow(regions), length(types), dimnames=list(rownames(regions), types))

for(x in 1:nrow(regions)){
  tbl <- table(mlist[[x]][, "gene_biotype"])
  mm[x, names(tbl)] <- tbl
}

nEntries <- apply(mm,2,sum)
mm <- mm[,names(sort(nEntries, decreasing = TRUE))]

write.table(mm, "FeaturesVitaLoci.txt",sep="\t", quote=FALSE)

nSum <- apply(mm,1,sum)

mmF <- cbind(HIGA = 0, HI = 0, MOGA = 0, MO = 0, SUM = nSum, mm)

for(x in 1:nrow(regions)){
  mVEP <- read.csv(paste0("",rownames(regions)[x], ".snps.vep"),sep="\t", skip=86, header=FALSE)
  mGenAge <- read.csv(paste0("genes/all/",rownames(regions)[x], "_genAge.txt"), sep="\t", header=TRUE, row.names = 1)
  colnames(mVEP) <- c("Uploaded_variation","Location","Allele","Gene","Feature","Feature_type","Consequence",
                      "cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","Extra")

  cat("loaded\n")

  isInteresting <- c(grep("HIGH", mVEP[,"Extra"]), grep("MODERATE", mVEP[,"Extra"]))
  mSub <- mVEP[isInteresting,]
  genes <- unique(mlist[[x]][which(mlist[[x]][, "gene_biotype"] == "protein_coding"),1])

  conversion <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"), 
                      filters = "ensembl_gene_id", values = unique(genes), mart = mart)
  mSub <- cbind(mSub, "ensembl_gene_id" = NA)
  for(y in 1:nrow(mSub)){
    splitt <- strsplit(mSub[y, "Extra"], ";")
    n <- grep("SYMBOL=", splitt[[1]])
    g <- gsub("SYMBOL=", "", splitt[[1]][n])
    iix <- which(conversion[,3] == g)
    if(length(iix) == 1){
      mSub[y, "ensembl_gene_id"] <- conversion[which(conversion[,3] == g), 2]
    }
  }
  cat("conversion added\n")

  mm <- matrix(NA, length(genes), 13, dimnames = list(genes, c("external_gene_name", "Impact", "inGenAge", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "description", "genAge", "nLocations", "nFeatures", "variations")))
  for(gene in genes){
    #cat("Doing gene",gene,"\n")
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
    #cat("Step1\n")
    mm[gene, "nLocations"] <- length(unique(mSub[which(mSub[, "ensembl_gene_id"] == gene),"Location"]))
    mm[gene, "nFeatures"] <- length(unique(mSub[which(mSub[, "ensembl_gene_id"] == gene),"Feature"]))
    mm[gene, "variations"] <- paste0(mSub[which(mSub[, "ensembl_gene_id"] == gene),"Consequence"], collapse=";")
    #cat("Step2\n")
    isHIGH <- any(grepl("HIGH", mSub[which(mSub[, "ensembl_gene_id"] == gene),"Extra"]))
    isMODERATE <- any(grepl("MODERATE", mSub[which(mSub[, "ensembl_gene_id"] == gene),"Extra"]))

    mm[gene, "Impact"] <- "NONE"
    if(isMODERATE) mm[gene, "Impact"] <- "MODERATE"
    if(isHIGH) mm[gene, "Impact"] <- "HIGH"

    #cat("Step3\n")
    if(gene %in% rownames(mGenAge)){
      mm[gene, "inGenAge"] <- "YES"
      mm[gene, "genAge"] <- mGenAge[gene, "genAge"]
    }else{
      mm[gene, "inGenAge"] <- "NO"
    }
    #cat("Step4\n")
  }
  cat("Summarizing\n")
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
  cat("Writing\n")
  write.table(cbind(ensembl_gene_id = rownames(mm), mm), file = paste0("summary/",rownames(regions)[x], ".summary.txt"), 
              sep = "\t", quote = FALSE,na="", row.names=FALSE)
  write.table(cbind(region = rownames(regions)[x], cbind(ensembl_gene_id = rownames(mm), mm)), file = paste0("summary/combined.summary.txt"), 
              sep = "\t", quote = FALSE,na="", row.names=FALSE, append = TRUE, col.names = FALSE)
  mmF[rownames(regions)[x], "MO"] <- length(which(mm[,"Impact"] != "NONE"))
  mmF[rownames(regions)[x], "HI"] <- length(which(mm[,"Impact"] == "HIGH"))
  mmF[rownames(regions)[x], "MOGA"] <- length(which(mm[,"inGenAge"] == "YES"))
  mmF[rownames(regions)[x], "HIGA"] <- length(which(mm[,"Impact"] == "HIGH" & mm[,"inGenAge"] == "YES"))
}
write.table(mmF, "FeaturesVitaLociImpacts.txt",sep="\t", quote=FALSE)


