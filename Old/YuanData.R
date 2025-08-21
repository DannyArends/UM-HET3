setwd("C:/Users/Danny/Downloads")

mdata <- read.table("animalvals.csv", sep = ",", header= TRUE)

means.fem <- c()
sd.fem <- c()
means.mal <- c()
sd.mal <- c()
for(x in unique(mdata[, "strain"])){
  means.fem <- c(means.fem, mean(mdata[which(mdata[, "strain"] == x & mdata[, "sex"] == "f"), "value"] ))
  sd.fem <- c(sd.fem, sd(mdata[which(mdata[, "strain"] == x & mdata[, "sex"] == "f"), "value"] ))
  means.mal <- c(means.mal, mean(mdata[which(mdata[, "strain"] == x & mdata[, "sex"] == "m"), "value"] ))
  sd.mal <- c(sd.mal, sd(mdata[which(mdata[, "strain"] == x & mdata[, "sex"] == "m"), "value"] ))
}

plot(means.fem, means.mal, pch = 18)

# Add error bars for males and females
# Label the 4 UM-HET3 parental strains


