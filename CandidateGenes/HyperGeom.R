setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/GeneOntology/Raw Data")

aa <-read.table("KEGG_WG_summary.txt", row.names = NULL, sep = "\t")

isVita <- grep("Vita", aa[,2])
isSoma <- grep("Soma", aa[,2])

# Define our foreground (Significant DEG) and background (All)
fg <- aa[isVita, "subcategory"]
bg <- aa[, "subcategory"]

## Over / Under representation analysis of Gene Ontology
ORA <- c()
for (go in unique(fg)) {
  x <- length(which(fg == go))                        # FG annotations
  m <- length(which(bg == go))                        # BG annoattions
  n <- length(bg) - m                                 # BG without this annotation
  k <- length(fg)                                    # All FG probes

  ORA <- rbind(ORA, c(x, k, m, n + m, phyper(x - 1, m, n, k, lower.tail=FALSE)))
}
rownames(ORA) <- unique(fg)


# Define our foreground (Significant DEG) and background (All)
fg <- aa[isSoma, "subcategory"]
bg <- aa[, "subcategory"]

## Over / Under representation analysis of Gene Ontology
ORA <- c()
for (go in unique(fg)) {
  x <- length(which(fg == go))                        # FG annotations
  m <- length(which(bg == go))                        # BG annoattions
  n <- length(bg) - m                                 # BG without this annotation
  k <- length(fg)                                    # All FG probes

  ORA <- rbind(ORA, c(x, k, m, n + m, phyper(x - 1, m, n, k, lower.tail=FALSE)))
}

rownames(ORA) <- unique(fg)



library(readxl)
setwd("/home/rqdt9/Dropbox (UTHSC GGI)/ITP_HET3_Mapping_Paper_Arends_2021/GeneOntology")
file_path <- "WholeGenome_GeneOntology.xlsx"

# Read GO sheets
go_bp <- read_excel(file_path, sheet = "GO - BP")
go_mf <- read_excel(file_path, sheet = "GO - MF")
go_cc <- read_excel(file_path, sheet = "GO - CC")
kegg <- read_excel(file_path, sheet = "KEGG")
reactome <- read_excel(file_path, sheet = "REACTOME")


# Label categories
go_bp$Category <- "BP"
go_mf$Category <- "MF"
go_cc$Category <- "CC"
kegg$Category <- "KEGG"
reactome$Category <- "REACTOME"

headers <- unique(c(colnames(go_bp), colnames(go_mf), colnames(go_cc), colnames(kegg), colnames(reactome)))
rows <- nrow(go_bp) + nrow(go_mf) + nrow(go_cc) + nrow(kegg) + nrow(reactome)
mdata <- matrix(NA, rows, length(headers))
colnames(mdata) <- headers
ii <- 1
for(m in list(go_bp, go_mf, go_cc, kegg, reactome)){
  for(x in 1:nrow(m)){
    mdata[ii, colnames(m)] <- as.character(m[x,])
    ii <- ii + 1
  }
}
mdata <- data.frame(mdata)
# Minimal 3 (or more in the FG)
i <- which(as.numeric(unlist(lapply(strsplit(mdata[, "GeneRatio"], "/"), "[", 1))) >= 3)
mdata <- mdata[i, ]

# Minimum 5 fold enrichted
mdata <- mdata[which(as.numeric(mdata[, "FoldEnrichment"]) > 5),]

# Only include Loci with at least 3 entries
mdata <- mdata[which(mdata[, "Locus"] %in% names(which(table(mdata[,"Locus"]) > 2))),]

# Sort by enrichment
mdata <- mdata[sort(as.numeric(mdata[,"FoldEnrichment"]), index.return = TRUE, decreasing = TRUE)$ix,]

All_Genes_FC <- setNames(as.numeric(mdata$FoldEnrichment), 1:nrow(mdata))
GO_Terms_List <- split(1:nrow(mdata), mdata$Description)

min_size <- 10
max_size <- 500
GO_Terms_List <- GO_Terms_List[sapply(GO_Terms_List, length) >= min_size & sapply(GO_Terms_List, length) <= max_size]

# 3. Create the GO_Terms_List
# Group the new Gene IDs by their GO term description
GO_Terms_List <- split(mdata$GeneID, mdata$Description)

library(org.Mm.eg.db)
library(GOSemSim)
library(rrvgo)
semData <- godata('org.Mm.eg.db', ont="BP", computeIC=TRUE) 

simMatrix <- mgoSim(
    GO1 = mdata[which(mdata[, "Category"] == "BP"), "Term.ID"],
    GO2 = mdata[which(mdata[, "Category"] == "BP"), "Term.ID"],
    semData = semData, 
    measure = "Wang",
    combine = NULL
)

reduced_terms <- reduceSimMatrix(simMatrix = simMatrix, threshold = 0.7, gene2GO = mdata[which(mdata[, "Category"] == "BP"), "Term.ID"])

