# Manually merge regions keeping the minimal interval

mregions <- read.table("regions_4way_merged.txt", sep="\t")
for(x in 1:nrow(mregions)){
  cat(mregions[x,2], "\t", mregions[x,3], "\t", mregions[x,4], file = paste0("regionsFiles/",mregions[x,1], ".txt"))
}

# Call SNPs

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern, wait = FALSE)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}

callSNPs <- function(bamfiles, region = "Vita1a.txt") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/naszb2/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -R regions/", region," -a FORMAT/DP -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -i '%QUAL>=30 && INFO/DP>10' - -o ~/SNPVEPFeb23/",region,".vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

setwd("/home/danny/SNPVEPFeb23/")
regions <- list.files("regions")

for(r in regions){
  callSNPs(c("/naszb2/Mouse/DNA/Sequencing/UM-HET3/BALBBYJ_M_S2_L002_default_GRCm38_68.sorted.bam",
             "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
             "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam",
             "/naszb2/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam"
          ), r)
}

#### VEP

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern, wait = FALSE)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}

setwd("/home/danny/SNPVEPFeb23/")
regions <- list.files("regions")

for(r in regions){
  cmd <- paste0("/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPVEPFeb23/",region,".vcf --everything -o /home/danny/SNPVEPFeb23/",region,".vep")
  execute(cmd)
}
