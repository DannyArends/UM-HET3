setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

library(biomaRt)

bio.mart <- useMart("ENSEMBL_MART_ENSEMBL", host="https://nov2020.archive.ensembl.org", dataset="mmusculus_gene_ensembl")

regions <- c( "1:0:24042124",                   # Vita1a
              "1:7969906:34173196",             # Vita1b
              "1:111859457:133622154",          # Vita1c

              "2:111134233:120612428",          # Vita2a
              "2:132808839:140168146",          # Vita2b
              "2:154229551:158377039",          # Vita2c

              "3:42232742:121463445",           # Vita3a
              "3:68084635:121463445",           # Vita3b

              "4:46120492:55012301",            # Vita4a
              "4:66839598:86583085",            # Vita4b
              "4:144622766:156860686",          # Vita4c

              "5:36613621:90217582",            # Vita5a

              "6:48656204:62673982",            # Vita6a
              "6:97252294:119943820",           # Vita6b

              "8:25964832:55830812",            # Vita8a

              "9:13442519:39046745",            # Vita9a
              "9:54904313:58081975",            # Vita9b
              "9:95323685:124029281",           # Vita9c
              "9:124029281:124359700",          # Vita9d

              "10:56745317:100488594",          #Vita10a

              "11:0:30104580",                  #Vita11a
              "11:59009193:109367495",          #Vita11b

              "12:99576264:112855820",          #Vita12a

              "13:4574227:22278588",            #Vita13a
              "13:76135291:106345461",          #Vita13b

              "14:78415875:118874224",          #Vita14a

              "15:55481391:84914043",           #Vita15a

              "17:0:26542857",                  #Vita17a
              "17:8288690:25122899",            #Vita17b
              "17:53679006:78682893",           #Vita17c

              "18:46508371:90720763",           #Vita18a

              "X:0:69830745",                   #VitaXa
              "X:150646933:169476592")          #VitaXb

names(regions) <- c("Vita1a","Vita1b","Vita1c","Vita2a","Vita2b","Vita2c","Vita3a","Vita3b","Vita4a","Vita4b","Vita4c","Vita5a","Vita6a","Vita6b",
                    "Vita8a","Vita9a","Vita9b","Vita9c","Vita9d","Vita10a","Vita11a","Vita11b","Vita12a","Vita13a","Vita13b","Vita14a","Vita15a",
                    "Vita17a","Vita17b","Vita17c","Vita18a","VitaXa","VitaXb")

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
  reference = "/naszb2/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ",region," -a FORMAT/DP -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -i '%QUAL>=30 && INFO/DP>10' - -o ~/March14/", output, ".snps.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

bamfiles <- c("/naszb2/Mouse/DNA/Sequencing/UM-HET3/BALBBYJ_M_S2_L002_default_GRCm38_68.sorted.bam",
              "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
              "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam",
              "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam")

for(r in names(regions)[c(1,7)]){
  rr <- unlist(strsplit(regions[r],":"))
  if(rr[2] == 0) rr[2] <- 1
  reg <- paste0(rr[1],":", rr[2], "-", rr[3])
  callSNPs(bamfiles, reg, r)
}





setwd("~/March14")

files <- list.files(pattern = "\\.vcf$")

for(file in files){
## Do VEP predictions for all Files from david
cat(paste0("/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/", file, " --everything -o /home/danny/March14/", gsub("vcf","vep", file), "\n"))
}

# Commands generated
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita10a.snps.vcf --everything -o /home/danny/March14/Vita10a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita11a.snps.vcf --everything -o /home/danny/March14/Vita11a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita11b.snps.vcf --everything -o /home/danny/March14/Vita11b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita12a.snps.vcf --everything -o /home/danny/March14/Vita12a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita13a.snps.vcf --everything -o /home/danny/March14/Vita13a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita13b.snps.vcf --everything -o /home/danny/March14/Vita13b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita14a.snps.vcf --everything -o /home/danny/March14/Vita14a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita15a.snps.vcf --everything -o /home/danny/March14/Vita15a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita17a.snps.vcf --everything -o /home/danny/March14/Vita17a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita17b.snps.vcf --everything -o /home/danny/March14/Vita17b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita17c.snps.vcf --everything -o /home/danny/March14/Vita17c.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita18a.snps.vcf --everything -o /home/danny/March14/Vita18a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita1b.snps.vcf --everything -o /home/danny/March14/Vita1b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita1c.snps.vcf --everything -o /home/danny/March14/Vita1c.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita2a.snps.vcf --everything -o /home/danny/March14/Vita2a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita2b.snps.vcf --everything -o /home/danny/March14/Vita2b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita2c.snps.vcf --everything -o /home/danny/March14/Vita2c.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita3b.snps.vcf --everything -o /home/danny/March14/Vita3b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita4a.snps.vcf --everything -o /home/danny/March14/Vita4a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita4b.snps.vcf --everything -o /home/danny/March14/Vita4b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita4c.snps.vcf --everything -o /home/danny/March14/Vita4c.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita5a.snps.vcf --everything -o /home/danny/March14/Vita5a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita6a.snps.vcf --everything -o /home/danny/March14/Vita6a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita6b.snps.vcf --everything -o /home/danny/March14/Vita6b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita8a.snps.vcf --everything -o /home/danny/March14/Vita8a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita9a.snps.vcf --everything -o /home/danny/March14/Vita9a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita9b.snps.vcf --everything -o /home/danny/March14/Vita9b.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita9c.snps.vcf --everything -o /home/danny/March14/Vita9c.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita9d.snps.vcf --everything -o /home/danny/March14/Vita9d.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/VitaXa.snps.vcf --everything -o /home/danny/March14/VitaXa.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/VitaXb.snps.vcf --everything -o /home/danny/March14/VitaXb.snps.vep &

/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita1a.snps.vcf --everything -o /home/danny/March14/Vita1a.snps.vep &
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/March14/Vita3a.snps.vcf --everything -o /home/danny/March14/Vita3a.snps.vep &


setwd("/home/rqdt9/Dropbox (UTHSC GGI)/MyFolder/UM-HET3/March14")
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


