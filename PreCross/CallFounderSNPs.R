#
# callFounderSNPs.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# Code that runs SNP calling on the 4 founders per chromosome and performs VEP predictions
#

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern, wait = FALSE)
    cat(">>>>", res[1], "\n")
  }
}

callSNPs <- function(bamfiles, region = "1") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "~/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ",region," -a FORMAT/DP -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -i 'QUAL>=30 && INFO/DP>10' - | gzip > ~/UMHET3/", region, ".snps.vcf.gz")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

bamfiles <- c("~/NAS/Mouse/DNA/Sequencing/UM-HET3/BALBBYJ_M_S2_L002_default_GRCm38_68.sorted.bam",
              "~/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
              "~/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam",
              "~/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam")

for(chr in c(1:19, "X", "Y", "MT")) {
  callSNPs(bamfiles, chr)
}

setwd("~/UMHET3/")
vep <- "/home/danny/Github/ensembl-vep/vep"

files <- list.files(pattern = "\\.vcf.gz$")
for(file in files){
  execute(paste0(vep, " --species mus_musculus --refseq --offline -i ~/UMHET3/", file, " --everything -o ~/UMHET3/", gsub("vcf.gz", "vep", file)))
}

