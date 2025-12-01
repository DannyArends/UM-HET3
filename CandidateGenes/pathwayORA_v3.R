#
# pathwayORA_v3.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Updated version of pathway ORA analysis, We remove GPCRs and Olfactory genes
# Performs a leave-one-Vita region out approach to see if in-Region clustering of genes explains the observed over-representation
#

library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)

# Manually merge regions keeping the minimal interval
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3")
regions <- read.table("regions_Aug25.txt", sep="\t", row.names=1)
colnames(regions) <- c("Chr", "Proximal", "Distal")
regions[,"Proximal"] <- 1e6 *as.numeric(regions[,"Proximal"])
regions[,"Distal"] <- 1e6 *as.numeric(regions[,"Distal"])

GPCRs <- read.table("May2024/GPCRs.txt", sep = "\t", header=TRUE)

### Get whole genome
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org")
mlist <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "gene_biotype", 
                                     "chromosome_name", "start_position", "end_position", "strand", "description", "mgi_description"), mart = mart)

#iix <- which(mlist[, "ensembl_gene_id"] %in% GPCRs[,"Ensembl.Gene.ID..Mouse."])
#mlist <- mlist[-iix,]
olf <- grep("^Olf", mlist[,"external_gene_name"])
if(length(olf)> 0) mlist <- mlist[-olf,]

vnm <- grep("^Vmn", mlist[,"external_gene_name"])
if(length(vnm)> 0) mlist <- mlist[-vnm,]

mlist <- mlist[which(mlist[, "chromosome_name"] %in% c(1:19, "X")),]
mlist <- mlist[which(mlist[, "gene_biotype"] == "protein_coding"),]
isGM <- grep("^Gm", mlist[,"external_gene_name"])
mlist <- mlist[-isGM, ]
isRIKEN <- grep("^RIKEN", mlist[,"mgi_description"])
genome <- mlist[-isRIKEN, ]

## Convert from ensembl ID to entrez ID
getEntrez <- function(ensembl) { return(as.character(na.omit(genome[which(genome[,"ensembl_gene_id"] %in% ensembl), "entrezgene_id"]))) }
wholeGenome <- unique(genome[,"ensembl_gene_id"])


setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/regionsAug25")
protein_coding <- vector("list", nrow(regions))
prioritized <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){
  protein_coding[[x]] <- read.csv(paste0("genes/protein_coding_genes_", rownames(regions)[x], ".txt"), sep = "\t")
  #iix <- which(protein_coding[[x]][, "ensembl_gene_id"] %in% GPCRs[,"Ensembl.Gene.ID..Mouse."])
  #protein_coding[[x]] <- protein_coding[[x]][-iix,]
  olf <- grep("^Olf", protein_coding[[x]][,"external_gene_name"])
  if(length(olf)> 0) protein_coding[[x]] <- protein_coding[[x]][-olf,]

  vnm <- grep("^Vmn", protein_coding[[x]][,"external_gene_name"])
  if(length(vnm)> 0) protein_coding[[x]] <- protein_coding[[x]][-vnm,]

  prioritized[[x]] <- read.csv(paste0("summary/", rownames(regions)[x], ".summary.txt"), sep = "\t")
  #iix <- which(prioritized[[x]][, "ensembl_gene_id"] %in% GPCRs[,"Ensembl.Gene.ID..Mouse."])
  #prioritized[[x]] <- prioritized[[x]][-iix,]
  olf <- grep("^Olf", prioritized[[x]][,"name"])
  if(length(olf)> 0) protein_coding[[x]] <- prioritized[[x]][-olf,]

  vnm <- grep("^Vmn", prioritized[[x]][,"name"])
  if(length(vnm)> 0) protein_coding[[x]] <- prioritized[[x]][-vnm,]
}

## One region OR versus 
enrichments.BP <- vector("list", nrow(regions))
enrichments.MF <- vector("list", nrow(regions))
enrichments.CC <- vector("list", nrow(regions))

KEGGs <- vector("list", nrow(regions))
REACTOMEs <- vector("list", nrow(regions))
BGs <- vector("list", nrow(regions))
FGs <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){
  BGs[[x]] <- wholeGenome
  FGs[[x]] <- prioritized[[x]][, 1]

  # Gene Ontology
  enrichments.BP[[x]] <- enrichGO(FGs[[x]], "org.Mm.eg.db", ont = "BP", keyType="ENSEMBL", pvalueCutoff = 0.2, universe = BGs[[x]])
  enrichments.MF[[x]] <- enrichGO(FGs[[x]], "org.Mm.eg.db", ont = "MF", keyType="ENSEMBL", pvalueCutoff = 0.2, universe = BGs[[x]])
  enrichments.CC[[x]] <- enrichGO(FGs[[x]], "org.Mm.eg.db", ont = "CC", keyType="ENSEMBL", pvalueCutoff = 0.2, universe = BGs[[x]])

  # Kegg
  KEGGs[[x]] <- enrichKEGG(gene = getEntrez(FGs[[x]]), organism = 'mmu', pvalueCutoff = 0.2, universe = getEntrez(BGs[[x]]))
  KEGGs[[x]] <- KEGGs[[x]][, 1:9]

  # Reactome
  REACTOMEs[[x]] <- enrichPathway(gene = getEntrez(FGs[[x]]), organism = "mouse", pvalueCutoff = 0.2, universe = getEntrez(BGs[[x]]))
  REACTOMEs[[x]] <- REACTOMEs[[x]][, 1:7]
  cat("Finished with",x,"\n")
}

names(enrichments.BP) <- rownames(regions)
names(enrichments.MF) <- rownames(regions)
names(enrichments.CC) <- rownames(regions)
names(KEGGs) <- rownames(regions)
names(REACTOMEs) <- rownames(regions)

for(x in rownames(regions)){
  write.table(enrichments.BP[[x]], file = paste0("GO_WG_BP_",x,".txt"), sep = "\t", quote = FALSE)
  write.table(enrichments.MF[[x]], file = paste0("GO_WG_MF_",x,".txt"), sep = "\t", quote = FALSE)
  write.table(enrichments.CC[[x]], file = paste0("GO_WG_CC_",x,".txt"), sep = "\t", quote = FALSE)

  write.table(KEGGs[[x]], file = paste0("KEGG_WG_",x,".txt"), sep = "\t", quote = FALSE)
  write.table(REACTOMEs[[x]], file = paste0("REACTOME_WG_",x,".txt"), sep = "\t", quote = FALSE)
}

header <- TRUE
for(x in rownames(regions)){
  if(is.null(enrichments.BP[[x]][,1:10])) next;
  if(nrow(enrichments.BP[[x]][,1:10]) > 0){
    write.table(cbind(x, enrichments.BP[[x]][,1:10]), file = paste0("GO_BP_WG_summary.txt"), sep = "\t", quote = FALSE, append = TRUE, col.names = header)
    header <- FALSE
  }
}


header <- TRUE
for(x in rownames(regions)){
  if(is.null(enrichments.MF[[x]][,1:10])) next;
  if(nrow(enrichments.MF[[x]][,1:10]) > 0){
    write.table(cbind(x, enrichments.MF[[x]][,1:10]), file = paste0("GO_MF_WG_summary.txt"), sep = "\t", quote = FALSE, append = TRUE, col.names = header)
    header <- FALSE
  }
}

header <- TRUE
for(x in rownames(regions)){
  if(is.null(enrichments.CC[[x]][,1:10])) next;
  if(nrow(enrichments.CC[[x]][,1:10]) > 0){
    write.table(cbind(x, enrichments.CC[[x]][,1:10]), file = paste0("GO_CC_WG_summary.txt"), sep = "\t", quote = FALSE, append = TRUE, col.names = header)
    header <- FALSE
  }
}

header <- TRUE
for(x in rownames(regions)){
  if(is.null(KEGGs[[x]])) next;
  if(nrow(KEGGs[[x]]) > 0){
    write.table(cbind(x, KEGGs[[x]]), file = paste0("KEGG_WG_summary.txt"), sep = "\t", quote = FALSE, append = TRUE, col.names = header)
    header <- FALSE
  }
}

header <- TRUE
for(x in rownames(regions)){
  if(is.null(REACTOMEs[[x]])) next;
  if(nrow(REACTOMEs[[x]]) > 0){
    write.table(cbind(x, REACTOMEs[[x]]), file = paste0("REACTOME_WG_summary.txt"), sep = "\t", quote = FALSE, append = TRUE, col.names = header)
    header <- FALSE
  }
}

