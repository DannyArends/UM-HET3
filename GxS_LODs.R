
setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

all <- c("1_3010272", "1_24042124", "1_120474787", "2_89844287", "2_112712327", "2_148442635","3_83838529", "3_92135706", "4_52524395",
         "5_67573068", "6_107382038", "6_132762500", "9_29939029", "9_104091597", "9_124056586", "10_72780332", "11_5628810", "11_82178599",
         "12_112855820", "13_89689878", "14_101437457", "15_74248242", "17_32883804", "18_60822951", "X_36008085", "X_156343080")

lods.i <- read.table("progressiveMapping_INT.txt", sep = "\t", check.names = FALSE)

lods.i[, all]
