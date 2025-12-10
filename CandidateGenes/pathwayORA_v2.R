#
# pathwayORA_v2.R
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

setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3")

regions <- read.table("regions_4way_topm.txt", sep="\t", header=FALSE, row.names=1)
colnames(regions) <- c("Chr", "Proximal", "Top", "Distal")
regions <- regions[which(!(rownames(regions) %in% c("Vita17a", "VitaXa", "VitaXb"))),]


GPCRs <- read.table("May2024/GPCRs.txt", sep = "\t", header=TRUE)

### Get whole genome
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org")
mlist <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "gene_biotype", 
                                     "chromosome_name", "start_position", "end_position", "strand", "description", "mgi_description"), mart = mart)

iix <- which(mlist[, "ensembl_gene_id"] %in% GPCRs[,"Ensembl.Gene.ID..Mouse."])
mlist <- mlist[-iix,]
olf <- grep("^Olf", mlist[,"external_gene_name"])
if(length(olf)> 0) mlist <- mlist[-olf,]
mlist <- mlist[which(mlist[, "chromosome_name"] %in% c(1:19, "X")),]
mlist <- mlist[which(mlist[, "gene_biotype"] == "protein_coding"),]
isGM <- grep("^Gm", mlist[,"external_gene_name"])
mlist <- mlist[-isGM, ]
isRIKEN <- grep("^RIKEN", mlist[,"mgi_description"])
genome <- mlist[-isRIKEN, ]

## Convert from ensembl ID to entrez ID

getEntrez <- function(ensembl) { return(as.character(na.omit(genome[which(genome[,"ensembl_gene_id"] %in% ensembl), "entrezgene_id"]))) }
wholeGenome <- unique(genome[,"ensembl_gene_id"])

## Load our data
protein_coding <- vector("list", nrow(regions))
prioritized <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){
  protein_coding[[x]] <- read.csv(paste0("2025/genes/protein_coding_genes_", rownames(regions)[x], ".txt"), sep = "\t")
  iix <- which(protein_coding[[x]][, "ensembl_gene_id"] %in% GPCRs[,"Ensembl.Gene.ID..Mouse."])
  protein_coding[[x]] <- protein_coding[[x]][-iix,]
  olf <- grep("^Olf", protein_coding[[x]][,"external_gene_name"])
  if(length(olf)> 0) protein_coding[[x]] <- protein_coding[[x]][-olf,]

  # Distance max 2 Mb from top
  d <- abs(((protein_coding[[x]][, "start_position"] + protein_coding[[x]][, "end_position"]) / 2) - regions[x, "Top"])
  protein_coding[[x]] <- protein_coding[[x]][which(d < 1000000),]

  prioritized[[x]] <- read.csv(paste0("2025/summary/", rownames(regions)[x], ".summary.txt"), sep = "\t")
  iix <- which(prioritized[[x]][, "ensembl_gene_id"] %in% GPCRs[,"Ensembl.Gene.ID..Mouse."])
  prioritized[[x]] <- prioritized[[x]][-iix,]
  olf <- grep("^Olf", prioritized[[x]][,"name"])
  if(length(olf)> 0) protein_coding[[x]] <- prioritized[[x]][-olf,]

  # Distance max 2 Mb from top
  d <- abs(((prioritized[[x]][, "start"] + prioritized[[x]][, "stop"]) / 2) - regions[x, "Top"])
  prioritized[[x]] <- prioritized[[x]][which(d < 1000000),]
}

## Leave one region out
enrichments <- vector("list", nrow(regions))
KEGGs <- vector("list", nrow(regions))
REACTOMEs <- vector("list", nrow(regions))
BGs <- vector("list", nrow(regions))
FGs <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){
  bg <- c()
  fg <- c()
  for(y in (1:nrow(regions))[-x]){
    iiMod <- which(prioritized[[y]][,3] == "HIGH")
    bg <- c(bg, protein_coding[[y]][,1])
    fg <- c(fg, prioritized[[y]][iiMod,1])
  }
  BGs[[x]] <- bg
  FGs[[x]] <- fg
  # Gene Ontology
  enrichments[[x]] <- enrichGO(FGs[[x]], "org.Mm.eg.db", ont = "BP", keyType="ENSEMBL", pvalueCutoff = 0.2, universe = BGs[[x]])
  enrichments[[x]] <- enrichments[[x]][, 1:8]

  # Kegg
  KEGGs[[x]] <- enrichKEGG(gene = getEntrez(FGs[[x]]), organism = 'mmu', pvalueCutoff = 0.2, universe = getEntrez(BGs[[x]]))
  KEGGs[[x]] <- KEGGs[[x]][, 1:9]

  # Reactome
  REACTOMEs[[x]] <- enrichPathway(gene = getEntrez(FGs[[x]]), organism = "mouse", pvalueCutoff = 0.2, universe = getEntrez(BGs[[x]]))
  REACTOMEs[[x]] <- REACTOMEs[[x]][, 1:7]
  cat("Finished with",x,"\n")
}

names(enrichments) <- rownames(regions)
names(KEGGs) <- rownames(regions)
names(REACTOMEs) <- rownames(regions)

any(table(unlist(lapply(enrichments,rownames))) == nrow(regions))
any(table(unlist(lapply(KEGGs,rownames))) == nrow(regions))
any(table(unlist(lapply(REACTOMEs,rownames))) == nrow(regions))


