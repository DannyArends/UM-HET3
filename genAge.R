#
# genAge
#
setwd("C:/Users/rqdt9/OneDrive - Northumbria University - Production Azure AD/Documents/HU-Berlin/UM-HET3/files")

GAmodel <- read.table("genage/genage_models.csv",sep=",", header=TRUE)
GAhuman <- read.table("genage/genage_human.csv",sep=",", header=TRUE)

regions <- read.table("regions_pmap.txt", sep="\t", header=TRUE)

library(biomaRt)

mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org")

for(x in 1:nrow(regions)){
  r <- paste0(regions[x, "Chr"],":", regions[x, "Proximal"], ":", regions[x, "Distal"])
  ortho <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_associated_gene_name", "celegans_homolog_associated_gene_name", "dmelanogaster_homolog_associated_gene_name", "scerevisiae_homolog_associated_gene_name"), filters = "chromosomal_region", values = r, mart = mart)
  descr <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "gene_biotype", "mgi_id", "mgi_symbol", "mgi_description"), filters = "chromosomal_region", values = r, mart = mart)

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
  write.table(descr, paste0("genes/all/", gsub(":", "_", r), ".txt"), sep = "\t", quote = FALSE, na = "")
  cat("Done",r, "in GenAge:", length(unique(inGenAge[,2])), "/", sum(descr[, "gene_biotype"] == "protein_coding"), "total Features:", nrow(descr), "\n")
}

#Done 1:0:24042124 in GenAge: 9 / 81 total Features: 319 
#Done 1:7969906:34173196 in GenAge: 8 / 87 total Features: 373 
#Done 2:92975818:153078512 in GenAge: 36 / 533 total Features: 1068 
#Done 4:30761996:107374161 in GenAge: 34 / 480 total Features: 1150 
#Done 6:83151312:122290952 in GenAge: 24 / 304 total Features: 768 
#Done 9:95323685:124029281 in GenAge: 20 / 356 total Features: 809 
#Done 10:56745317:106030562 in GenAge: 48 / 496 total Features: 1151 
#Done 12:99576264:120092757 in GenAge: 12 / 180 total Features: 801 
#Done 13:76135291:106345461 in GenAge: 15 / 141 total Features: 430 
#Done 14:78415875:118874224 in GenAge: 4 / 62 total Features: 366 
#Done 15:55481391:88420584 in GenAge: 30 / 361 total Features: 765 
#Done 17:0:48308023 in GenAge: 61 / 809 total Features: 1609 
#Done 18:34938676:90720763 in GenAge: 28 / 377 total Features: 986 
#Done X:0:69830745 in GenAge: 22 / 376 total Features: 997 
#Done X:150646933:169476592 in GenAge: 11 / 114 total Features: 306 
#Done 1:7969906:34173196 in GenAge: 8 / 87 total Features: 373 
#Done 2:122685194:182113224 in GenAge: 45 / 674 total Features: 1294 
#Done 3:68084635:121463445 in GenAge: 38 / 569 total Features: 1295 
#Done 9:13442519:39046745 in GenAge: 14 / 291 total Features: 628 
#Done 9:58081975:124359700 in GenAge: 45 / 611 total Features: 1558 
#Done 11:59009193:109367495 in GenAge: 92 / 953 total Features: 1658 
#Done X:150646933:169476592 in GenAge: 11 / 114 total Features: 306 
#Done 1:0:24042124 in GenAge: 9 / 81 total Features: 319 
#Done 1:87986124:143770506 in GenAge: 14 / 305 total Features: 858 
#Done 2:71846197:139956785 in GenAge: 44 / 782 total Features: 1458 
#Done 2:71846197:139956785 in GenAge: 44 / 782 total Features: 1458 
#Done 3:42232742:121463445 in GenAge: 43 / 672 total Features: 1697 
#Done 4:30761996:63272545 in GenAge: 16 / 281 total Features: 623 
#Done 5:36613621:90217582 in GenAge: 14 / 232 total Features: 843 
#Done 6:97252294:149588044 in GenAge: 34 / 452 total Features: 1153 
#Done 9:108054895:124359700 in GenAge: 10 / 210 total Features: 461 
#Done 10:46831195:100488594 in GenAge: 50 / 504 total Features: 1225 
#Done 11:0:30104580 in GenAge: 17 / 165 total Features: 431 
#Done 15:34609241:84914043 in GenAge: 33 / 409 total Features: 969 
#Done 17:8288690:53515555 in GenAge: 59 / 792 total Features: 1580 

regions <- read.table("regions_pmap_ctrl.txt", sep="\t", header=TRUE)

library(biomaRt)

mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org")

for(x in 1:nrow(regions)){
  r <- paste0(regions[x, "Chr"],":", regions[x, "Proximal"], ":", regions[x, "Distal"])
  ortho <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_associated_gene_name", "celegans_homolog_associated_gene_name", "dmelanogaster_homolog_associated_gene_name", "scerevisiae_homolog_associated_gene_name"), filters = "chromosomal_region", values = r, mart = mart)
  descr <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "gene_biotype", "mgi_id", "mgi_symbol", "mgi_description"), filters = "chromosomal_region", values = r, mart = mart)

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
  write.table(descr, paste0("genes/ctrl/", gsub(":", "_", r), ".txt"), sep = "\t", quote = FALSE, na = "")
  cat("Done",r, "in GenAge:", length(unique(inGenAge[,2])), "/", sum(descr[, "gene_biotype"] == "protein_coding"), "total Features:", nrow(descr), "\n")
}

#Done 1:42044440:58823920 in GenAge: 18 / 87 total Features: 276 
#Done 2:30426723:129010510 in GenAge: 64 / 949 total Features: 1893 
#Done 3:89167102:121463445 in GenAge: 30 / 435 total Features: 950 
#Done 9:0:124359700 in GenAge: 86 / 1205 total Features: 2838 
#Done 9:58081975:124359700 in GenAge: 45 / 611 total Features: 1558 
#Done 12:99576264:120092757 in GenAge: 12 / 180 total Features: 801 
#Done 13:83858506:120883175 in GenAge: 20 / 200 total Features: 582 
#Done 14:48473970:106980923 in GenAge: 30 / 435 total Features: 1282 
#Done 17:0:32883804 in GenAge: 36 / 437 total Features: 952 
#Done 19:8803985:38057964 in GenAge: 22 / 297 total Features: 618 
#Done 2:92975818:182113224 in GenAge: 63 / 930 total Features: 1805 
#Done 3:93148504:121463445 in GenAge: 25 / 316 total Features: 768 
#Done 9:13442519:54904313 in GenAge: 32 / 498 total Features: 1076 
#Done 16:30343647:74899626 in GenAge: 15 / 248 total Features: 610 
#Done 19:4907657:43849879 in GenAge: 39 / 508 total Features: 976 
#Done 1:170283138:195154279 in GenAge: 18 / 244 total Features: 575 
#Done 2:26373532:119210795 in GenAge: 61 / 892 total Features: 1810 
#Done 3:28380854:121463445 in GenAge: 52 / 748 total Features: 1961 
#Done 4:134348530:156860686 in GenAge: 25 / 354 total Features: 669 
#Done 9:120427595:124359700 in GenAge: 4 / 56 total Features: 146 
#Done 10:56745317:122981479 in GenAge: 58 / 589 total Features: 1459 
#Done 12:79076511:114630804 in GenAge: 17 / 315 total Features: 1072 
#Done 14:0:41200458 in GenAge: 21 / 321 total Features: 769 
#Done 15:48895974:102208446 in GenAge: 49 / 589 total Features: 1208 
#Done 17:0:28902980 in GenAge: 29 / 374 total Features: 789 

