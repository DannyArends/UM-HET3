setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files/2025_Genes")

library(biomaRt)

bio.mart <- useMart("ENSEMBL_MART_ENSEMBL", host="https://nov2020.archive.ensembl.org", dataset="mmusculus_gene_ensembl")

regions <- c( "1:0:7969906",                    # Vita1a n
              "1:7969906:34173196",             # Vita1b n
              "1:87986124:143770506",           # Vita1c n

              "2:71846197:99777113",            # Vita2a n
              "2:92975818:125083168",           # Vita2b n
              "2:132808839:161871392",          # Vita2c n

              "3:42232742:121463445",           # Vita3a n
              "3:68084635:121463445",           # Vita3b n

              "4:30761996:63272545",            # Vita4a n

              "5:36613621:90217582",            # Vita5a n

              "6:97252294:122290952",           # Vita6a n
              "6:108075853:149736546",          # Vita6b n

              "9:13442519:44688426",            # Vita9a n
              "9:95323685:116536254",           # Vita9b n
              "9:108054895:124359700",          # Vita9c n

              "10:56745317:100488594",          #Vita10a n

              "11:0:30104580",                  #Vita11a n
              "11:59009193:109367495",          #Vita11b n

              "12:99576264:120129022",          #Vita12a n

              "13:76135291:104111134",          #Vita13a n

              "14:78415875:118874224",          #Vita14a n

              "15:62405371:88420584",           #Vita15a n

              "17:8288690:48308023",            #Vita17a n

              "18:38240954:70827226",           #Vita18a n

              "X:0:69830745",                   #VitaXa n
              "X:150646933:169476592")          #VitaXb n

names(regions) <- c("Vita1a","Vita1b","Vita1c","Vita2a","Vita2b","Vita2c","Vita3a","Vita3b","Vita4a","Vita5a","Vita6a","Vita6b",
                    "Vita9a","Vita9b","Vita9c","Vita10a","Vita11a","Vita11b","Vita12a","Vita13a","Vita14a","Vita15a",
                    "Vita17a","Vita18a","VitaXa","VitaXb")


res <- vector("list", length(regions))
cnt <- 1
for(r in regions){
  res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", 
                                      "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description", "gene_biotype"), 
                       filters = c("chromosomal_region"), values = list(r), mart = bio.mart)
  isGM <- grep("^Gm", res.biomart[,"external_gene_name"])
  if(length(isGM) > 0) res.biomart <- res.biomart[-isGM,]

  isRIKEN <- grep("^RIKEN", res.biomart[,"mgi_description"])
  if(length(isRIKEN) > 0) res.biomart <- res.biomart[-isRIKEN,]

  write.table(res.biomart, file=paste0(names(regions)[cnt], ".features.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  res[[cnt]] <- res.biomart
  cat("Done",r,"\n")
  cnt <- cnt + 1
}

nTotal <- unlist(lapply(lapply(res, dim),"[",1))
nGm <- unlist(lapply(lapply(res, function(x){length(grep("^Gm", x[,"external_gene_name"]))}),"[",1))
mTab <- lapply(res,function(x){
  table(x[,"gene_biotype"])
})
nBioTypes <- unique(unlist(lapply(mTab, names)))
mmatrix <- matrix(NA, length(regions), length(nBioTypes), dimnames = list(regions, nBioTypes))

for(x in 1:length(mTab)){
  mmatrix[x, names(mTab[[x]])] <- mTab[[x]]
}
mmatrix[is.na(mmatrix)] <- 0

nEntries <- apply(mmatrix,2,sum)

mmatrix <- mmatrix[,names(sort(nEntries, decreasing = TRUE))]
rownames(mmatrix) <- names(regions)

write.table(mmatrix, file = "GenBiotypes.txt",sep="\t",quote=FALSE,row.names=TRUE, na = "")

#### Berlin Server part

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern, wait = FALSE)
    cat(">>>>", res[1], "\n")
  }
}

callSNPs <- function(bamfiles, region = "chr1:124548738-184926264", output = "") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "~/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ",region," -a FORMAT/DP -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -i 'QUAL>=30 && INFO/DP>10' - -o ~/2025/", output, ".snps.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

bamfiles <- c("~/NAS/Mouse/DNA/Sequencing/UM-HET3/BALBBYJ_M_S2_L002_default_GRCm38_68.sorted.bam",
              "~/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
              "~/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam",
              "~/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam")

for(r in names(regions)){
  rr <- unlist(strsplit(regions[r],":"))
  if(rr[2] == 0) rr[2] <- 1
  reg <- paste0(rr[1],":", rr[2], "-", rr[3])
  callSNPs(bamfiles, reg, r)
}





setwd("~/2025")

files <- list.files(pattern = "\\.vcf$")

for(file in files){
## Do VEP predictions for all vcf files generated
cat(paste0("/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/", file, " --everything -o /home/danny/2025/", gsub("vcf","vep", file), "\n"))
}

# Commands generated
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita10a.snps.vcf --everything -o /home/danny/2025/Vita10a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita11a.snps.vcf --everything -o /home/danny/2025/Vita11a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita11b.snps.vcf --everything -o /home/danny/2025/Vita11b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita12a.snps.vcf --everything -o /home/danny/2025/Vita12a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita13a.snps.vcf --everything -o /home/danny/2025/Vita13a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita14a.snps.vcf --everything -o /home/danny/2025/Vita14a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita15a.snps.vcf --everything -o /home/danny/2025/Vita15a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita17a.snps.vcf --everything -o /home/danny/2025/Vita17a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita18a.snps.vcf --everything -o /home/danny/2025/Vita18a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita1a.snps.vcf --everything -o /home/danny/2025/Vita1a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita1b.snps.vcf --everything -o /home/danny/2025/Vita1b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita1c.snps.vcf --everything -o /home/danny/2025/Vita1c.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita2a.snps.vcf --everything -o /home/danny/2025/Vita2a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita2b.snps.vcf --everything -o /home/danny/2025/Vita2b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita2c.snps.vcf --everything -o /home/danny/2025/Vita2c.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita3a.snps.vcf --everything -o /home/danny/2025/Vita3a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita3b.snps.vcf --everything -o /home/danny/2025/Vita3b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita4a.snps.vcf --everything -o /home/danny/2025/Vita4a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita5a.snps.vcf --everything -o /home/danny/2025/Vita5a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita6a.snps.vcf --everything -o /home/danny/2025/Vita6a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita6b.snps.vcf --everything -o /home/danny/2025/Vita6b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita9a.snps.vcf --everything -o /home/danny/2025/Vita9a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita9b.snps.vcf --everything -o /home/danny/2025/Vita9b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/Vita9c.snps.vcf --everything -o /home/danny/2025/Vita9c.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/VitaXa.snps.vcf --everything -o /home/danny/2025/VitaXa.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --refseq --offline -i /home/danny/2025/VitaXb.snps.vcf --everything -o /home/danny/2025/VitaXb.snps.vep &




setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/2025")
files <- list.files(pattern = "\\.vep$")

mvepdata <- c()
for(f in files){
  mvepdata <- rbind(mvepdata, read.table(f))
}

mvepdata <- unique(mvepdata)
dim(mvepdata)

mvepdata <- mvepdata[grep("missense_variant", mvepdata[,"V7"]),]
dim(mvepdata)

for(r in regions){
  chr <- strsplit(r, ":")[[1]][1]; sta <- strsplit(r, ":")[[1]][2]; sto <- strsplit(r, ":")[[1]][3]
  ms <- strsplit(gsub("chr", "", mvepdata[,2]), ":")
  chrs <- unlist(lapply(ms, "[", 1))
  poss <- unlist(lapply(ms, "[", 2))
  iii <- which(chrs == chr & as.numeric(poss) > as.numeric(sta) & as.numeric(poss) < as.numeric(sto))
  nDel <- grep("SIFT=deleterious", mvepdata[iii, "V14"])

  namez <- mvepdata[iii, "V2"]
  namezDel <- mvepdata[iii[nDel], "V2"]
  namezGDel <- mvepdata[iii[nDel], "V4"]
  cat(r," ", length(unique(namez))," ", length(unique(namezDel))," ", length(unique(namezGDel)), "\n")
  if(length(namezGDel) > 0){
    res.biomart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "mgi_description"), 
                          filters = c("ensembl_gene_id"), values = list(unique(namezGDel)), mart = bio.mart)

    isGM <- grep("^Gm", res.biomart[,"external_gene_name"])
    if(length(isGM) > 0) res.biomart <- res.biomart[-isGM,]

    isRIKEN <- grep("^RIKEN", res.biomart[,"mgi_description"])
    if(length(isRIKEN) > 0) res.biomart <- res.biomart[-isRIKEN,]

    cat(paste0(res.biomart[, "external_gene_name"], collapse=", "), "\n")
  }
}


