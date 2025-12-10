install.packages("qtl", repos = "https://cloud.r-project.org")
install.packages("quantreg", repos = "https://cloud.r-project.org")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomaRt")

