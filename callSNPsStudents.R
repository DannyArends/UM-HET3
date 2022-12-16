callSNPs <- function(bamfiles, region = "chr1:124548738-184926264", output = "") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ",region," -a FORMAT/DP -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -i '%QUAL>=30 && INFO/DP>10' - -o ~/UMHET3parental/Founders.snps",output,".vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

bamfiles <- c("/home/danny/NAS/Mouse/DNA/Sequencing/UM-HET3/BALBBYJ_M_S2_L002_default_GRCm38_68.sorted.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam")

callSNPs(bamfiles, "chr18:11017447:36665388", "muqadassa")


/home/danny/Github/ensembl-vep/vep --species mus_musculus --offline -i /home/danny/SNPsDashbrook/Joint_called_gvcfs_ITP_parents_QTL_chr1_0_26682184.vcf  --everything -o /home/danny/SNPsDashbrook/Joint_called_gveps_ITP_parents_QTL_chr1_0_26682184.vep