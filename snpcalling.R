# Call SNPs on HET3

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern, wait = TRUE)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}

callSNPs <- function(chr = "chr1", outname = "snps") {
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/danny/UCSC_mm10.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup --threads 4 -q 30 -b bamfiles.txt -a FORMAT/DP -r ", chr, " -Ou -f ", reference)
  cmd2 <- paste0(bcftools, " call --threads 4 -mv -Ov ")
  cmd3 <- paste0(bcftools, " view --threads 4 -i '%QUAL>=100 && INFO/DP>10' - -o ~/", outname, "snps-filtered_",chr,".vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

setwd("/home/danny/NAS/Mouse/UMHET3bams")
bamfiles <- dir(".", "*.bam$")
cat(paste0("./", bamfiles, collapse="\n"), file="bamfiles.txt")

for(x in c("X")){
  callSNPs(paste0("chr", x), "UM_HET3Juli/")
}

# Call SNPs in the Founders based on the genome positions we found in the UM-HET3
execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}

callSNPs <- function(bamfiles) {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -R regionsITP.txt -a FORMAT/DP -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -i '%QUAL>=30 && INFO/DP>10' - -o ~/UMHET3parental/Founders.snps-filteredDP.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

callSNPs(c("/home/danny/NAS/Mouse/DNA/Sequencing/UM-HET3/BALBBYJ_M_S2_L002_default_GRCm38_68.sorted.bam",
           "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
           "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam",
           "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam"
          ))
