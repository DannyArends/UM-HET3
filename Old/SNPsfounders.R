# Call SNPs on UM-HET3 Founders
cat("/home/danny/NAS/Mouse/DNA/Sequencing/UM-HET3/BALBBYJ_M_S2_L002_default_GRCm38_68.sorted.bam\n", file = "bamfiles.txt")
cat("/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam\n", file = "bamfiles.txt", append = TRUE)
cat("/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam\n", file = "bamfiles.txt", append = TRUE)
cat("/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam", file = "bamfiles.txt", append = TRUE)

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern, wait = FALSE)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}


callSNPs <- function(chr = "chr1", outname = "snps") {
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup --threads 4 -q 30 -b bamfiles.txt -a FORMAT/DP -r ", chr, " -Ou -f ", reference)
  cmd2 <- paste0(bcftools, " call --threads 4 -mv -Ov ")
  cmd3 <- paste0(bcftools, " view --threads 4 -i 'QUAL>=100 && INFO/DP>10' - -o ", outname, "snps-filtered_",chr,".vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

for(x in c(1:19, "X")) {
  callSNPs(paste0(x), "~/wholeGenome/")
}


for(x in c(1:19, "X")) {
  cmd <- paste0("/home/danny/Github/ensembl-vep/vep --species mus_musculus_refseq --offline -i /home/danny/wholeGenome/snps-filtered_",x,".vcf --everything -o /home/danny/wholeGenome/VEP_chr",x,".vep")
  execute(cmd)
}
