setwd("C:/Github/UM-HET3/files")

library(biomaRt)

bio.mart <- useMart("ENSEMBL_MART_ENSEMBL", host="https://nov2020.archive.ensembl.org", dataset="mmusculus_gene_ensembl")

regions <- c( "1:0:26682184", 
              "1:7969906:31062737", 
              "2:122685194:153078512",
              "4:43036567:66839598",
              "4:74811205:112765106",
              "6:97252294:122290952",
              "9:13442519:39046745",
              "9:101942923:124359700",
              "9:90405149:124359700",
              "10:56745317:94163743",
              "12:99576264:120092757",
              "14:78415875:120296832",
              "15:48895974:84914043",
              "17:8288690:48308023",
              "X:0:69830745",
              "X:150646933:169476592")

res <- vector("list", length(regions))
cnt <- 1
for(r in regions){
  res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description", "gene_biotype"), 
                          filters = c("chromosomal_region"), values = list(r), mart = bio.mart)
  write.table(res.biomart, file=paste0("Tllq_",cnt, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
  res[[cnt]] <- res.biomart
  cat("Done",r,"\n")
  cnt <- cnt + 1
}

nTotal <- unlist(lapply(lapply(res, dim),"[",1))
nGm <- unlist(lapply(lapply(res, function(x){length(grep("^Gm", x[,"external_gene_name"]))}),"[",1))
mTab <- lapply(res,function(x){
  isGM <- grep("^Gm", x[,"external_gene_name"])
  table(x[-isGM,"gene_biotype"])
})
nBioTypes <- unique(unlist(lapply(mTab, names)))
mmatrix <- matrix(NA, length(regions), length(nBioTypes), dimnames = list(regions, nBioTypes))

for(x in 1:length(mTab)){
  mmatrix[x, names(mTab[[x]])] <- mTab[[x]]
}
write.table(mmatrix, file = "GenBiotypes.txt",sep="\t",quote=FALSE,row.names=FALSE, na = "")



setwd("~/SNPsDashbrook")

files <- list.files(pattern = "\\.vcf$")

for(file in files){
## Do VEP predictions for all Files from david
cat("/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/", file, " --everything -o /home/danny/SNPsDashbrook/", gsub("vcf","vep", file), "\n")
}

# Commands generated
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr1_0_26682184.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr1_0_26682184.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr1_0_31062759.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr1_0_31062759.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr1_0_34173196.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr1_0_34173196.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr10_46831195_94163743.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr10_46831195_94163743.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr10_56745317_94163743.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr10_56745317_94163743.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr12_99576264_120092757.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr12_99576264_120092757.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr14_78415875_120296832.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr14_78415875_120296832.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr15_31335181_84914043.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr15_31335181_84914043.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr15_48895974_84914043.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr15_48895974_84914043.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr17_0_48308023.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr17_0_48308023.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr17_8288690_53515555.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr17_8288690_53515555.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr1_7969906_31062737.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr1_7969906_31062737.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr2_113364081_182113224.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr2_113364081_182113224.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr2_122685194_182113224.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr2_122685194_182113224.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr2_92975818_153078512.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr2_92975818_153078512.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr4_30761996_66839598.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr4_30761996_66839598.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr4_43036567_74811205.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr4_43036567_74811205.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr4_74811205_112765106.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr4_74811205_112765106.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr6_81943533_122290952.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr6_81943533_122290952.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr6_97252294_149736546.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr6_97252294_149736546.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr9_101942923_124359700.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr9_101942923_124359700.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr9_13442519_39046745.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr9_13442519_39046745.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr9_13442519_44688426.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr9_13442519_44688426.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr9_90405149_124359700.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr9_90405149_124359700.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chrX_0_69830745.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chrX_0_69830745.vep
/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chrX_150646933_169476592.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chrX_150646933_169476592.vep


setwd("C:/Github/UM-HET3/files/DGA_vcfs")
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
  res.biomart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                          filters = c("ensembl_gene_id"), values = list(unique(namezGDel)), mart = bio.mart)
  cat(paste0(res.biomart[, "external_gene_name"], collapse=", "), "\n")
}
