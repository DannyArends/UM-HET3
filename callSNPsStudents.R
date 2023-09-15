execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern, wait = TRUE)
    cat(">>>>", res[1], "\n")
  }
}

callSNPs <- function(bamfiles, region = "chr1:124548738-184926264", output = "") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/naszb2/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ",region," -a FORMAT/DP -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -i '%QUAL>=30 && INFO/DP>10' - -o ~/Solomon/Founders.snps.",output,".vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

bamfiles <- c("/naszb2/Mouse/DNA/Sequencing/UM-HET3/BALBBYJ_M_S2_L002_default_GRCm38_68.sorted.bam",
              "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
              "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam",
              "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam")

aaa <- read.table("Solomon_Regions.txt")

for(x in 1:nrow(aaa)){
  region <- paste0(aaa[x, 1], ":",aaa[x, 2], "-",aaa[x, 3])
  callSNPs(bamfiles, region, region)
  cmd <- paste0("/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/Solomon/Founders.snps.",
                region,".vcf --everything -o /home/danny/Solomon/Founders.snps.",region,".vep")
  execute(cmd)
}



/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/StudentsNCL/Founders.snps.muqadassa.vcf  --everything -o /home/danny/StudentsNCL/Founders.snps.muqadassa.vep

/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/StudentsNCL/Founders.snps.taya.vcf  --everything -o /home/danny/StudentsNCL/Founders.snps.taya.vep

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern, wait = TRUE)
    cat(">>>>", res[1], "\n")
  }
}

callSNPs <- function(bamfiles, region = "chr1:124548738-184926264", output = "") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/naszb2/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ",region," -a FORMAT/DP -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -i '%QUAL>=30 && INFO/DP>10' - -o ~/Oghenetejiri/Founders.snps.",output,".vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

bamfiles <- c("/naszb2/Mouse/DNA/Sequencing/BFMI860-S12_2019/BAM/merged_sorted_860-S12.bam",
              "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam")

aaa <- read.table("Regions_Oghenetejiri.txt")

for(x in 1:nrow(aaa)){
  region <- paste0(aaa[x, 1], ":",aaa[x, 2], "-",aaa[x, 3])
  callSNPs(bamfiles, region, region)
  cmd <- paste0("/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/Oghenetejiri/Founders.snps.",
                region,".vcf --everything -o /home/danny/Oghenetejiri/Founders.snps.",region,".vep")
  execute(cmd)
}


