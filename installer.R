install.packages("qtl", repos = "https://cloud.r-project.org")
install.packages("quantreg", repos = "https://cloud.r-project.org")
install.packages("RColorBrewer", repos = "https://cloud.r-project.org")
install.packages("svglite", repos = "https://cloud.r-project.org")
install.packages("pspline", repos = "https://cloud.r-project.org")
install.packages("vioplot", repos = "https://cloud.r-project.org")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomaRt")

