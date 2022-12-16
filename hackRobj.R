setwd("C:/Users/rqdt9/OneDrive - Northumbria University - Production Azure AD/Documents/HU-Berlin/UM-HET3/files")

csvr <- read.table("um-het3-rqtl.csvr",sep=',')

wF <- rbind(
  csvr[1:10,],
  c("1_0", "1", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 1), ],
  c("2_0", "2", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 2), ],
  c("3_0", "3", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 3), ],
  c("4_0", "4", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 4), ],
  c("5_0", "5", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 5), ],
  c("6_0", "6", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 6), ],
  c("7_0", "7", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 7), ],
  c("8_0", "8", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 8), ],
  c("9_0", "9", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 9), ],
  c("10_0", "10", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 10), ],
  c("11_0", "11", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 11), ],
  c("12_0", "12", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 12), ],
  c("13_0", "13", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 13), ],
  c("14_0", "14", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 14), ],
  c("15_0", "15", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 15), ],
  c("16_0", "16", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 16), ],
  c("17_0", "17", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 17), ],
  c("18_0", "18", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 18), ],
  c("19_0", "19", "0", rep("", 6438)),
  csvr[which(csvr[,2] == 19), ],
  c("X_0", "X", "0", rep("", 6438)),
  csvr[which(csvr[,2] == "X"), ])

write.table(wF, "um-het3-rqtl.fake.csvr",sep=',', quote=FALSE, row.names=FALSE, col.names=FALSE, na = "")

# Read cross object
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.fake.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 1)
mcross <- adjustXprobs(mcross)

aa <- scanone(mcross)
write.table(aa, "map.all.txt", sep = "\t", quote=FALSE)
