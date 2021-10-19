setwd("C:/Github/UM-HET3/files")

gts4way <- read.table("merged/gts4way.rqtl.GN2.Juli2021.txt", sep="\t")
map <- read.table("merged/map.gts4way.Juli2021.txt", sep="\t")
ind <- read.table("merged/ind.gts4way.Juli2021.txt", sep="\t")

#map <- data.frame(cbind(Chr = gts4way[, "Chr"], Position = gts4way[, "Position"] * 1000000))
gts4wayRqtl <- gts4way[, -c(1:2)]

gts4wayNum <- apply(gts4wayRqtl, 2, as.numeric)
rownames(gts4wayNum) <- rownames(gts4wayNum)
table(gts4wayNum)

#gts4wayNum.cor <- cor(t(gts4wayNum), use="pair")
#write.table(gts4wayNum.cor, "merged/gts4way.raw.cor.txt", quote=FALSE, sep="\t")

# Create cross object for determining genotypes
setwd("C:/Github/UM-HET3/files")
write.table(cbind(NA,NA, t(cbind(Individual = rownames(ind), ind))), file = "um-het3-rqtl.csvr", col.names = FALSE, sep = ",", quote=FALSE, na = "")
write.table(rbind(GenoID = c(NA, NA, colnames(gts4wayRqtl)), 
                  cbind(Chr = map$Chr, Mb = as.numeric(map$Position) / 1000000, gts4wayRqtl)
                 ), file = "um-het3-rqtl.csvr", col.names = FALSE, sep = ",", append=TRUE, quote=FALSE, na="")

# Read cross object for determining genotypes
library(qtl)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross)

# Sequenom data
sequ <- subset(mcross, ind = mcross$pheno$Origin == "Sequenom")
sequ <- calc.genoprob(sequ)

# Monsterplex data
mons <- subset(mcross, ind = mcross$pheno$Origin == "Monsterplex")
mons <- calc.genoprob(mons)

mcross <- calc.genoprob(mcross)
plot.map(mcross, main = "Physical positions", ylab="Location (Mb)")

fcross <- fill.geno(mcross, method = "maxmarginal", error.prob = 0.01, min.prob = 0.85)
gts4way.full <- t(pull.geno(fcross))

mcor <- cor(t(gts4way.full), use="pair")
image(1:nrow(mcor), 1:ncol(mcor), mcor)

colnames(gts4way.full) <- as.character(fcross$pheno$Individual)

write.table(gts4way.full, "merged/gts4way.rqtl.filled.Juli2021.txt", quote=FALSE, sep="\t", na="")

gts4way.fullGN2 <- cbind(map[rownames(gts4way.full),c(2,3)], gts4way.full)

write.table(gts4way.fullGN2, "merged/gts4way.filled.GN2.txt", quote=FALSE, sep="\t", na="")

gts4way.fullnum <- apply(gts4way.full, 2, as.numeric)
rownames(gts4way.fullnum) <- rownames(gts4way.fullnum)

gts4way.cor <- cor(t(gts4way.fullnum), use="pair")
write.table(gts4way.cor, "gts4way.filled.cor.txt", quote=FALSE, sep="\t")

Treatment <- mcross$pheno$Treatment_Effect
Site <- mcross$pheno$Site
Sex <- mcross$pheno$Sex

females.JL <- subset(mcross, ind = Sex == "F" & Site == "JL" & Treatment %in% c("Cont", "Controls", "Aspirin"))
males.JL <- subset(mcross, ind = Sex == "M" & Site == "JL" & Treatment %in% c("Cont", "Controls", "Aspirin"))
females.UT <- subset(mcross, ind = Sex == "F" & Site == "UT" & Treatment %in% c("Cont", "Controls", "Aspirin"))
males.UT <- subset(mcross, ind = Sex == "M" & Site == "UT" & Treatment %in% c("Cont", "Controls", "Aspirin"))
females.UM <- subset(mcross, ind = Sex == "F" & Site == "UM" & Treatment %in% c("Cont", "Controls", "Aspirin"))
males.UM <- subset(mcross, ind = Sex == "M" & Site == "UM" & Treatment %in% c("Cont", "Controls", "Aspirin")) 

phe <- "BodyWeight_HET3_ITP_6m"#"Longevity_HET3_ITP"

qtl.JL.f <- scanone(females.JL, pheno.col = phe, method="hk")
qtl.JL.m <- scanone(males.JL, pheno.col = phe, method="hk")
plot(qtl.JL.f, qtl.JL.m)

qtl.UT.f <- scanone(females.UT, pheno.col = phe, method="hk")
qtl.UT.m <- scanone(males.UT, pheno.col = phe, method="hk")
plot(qtl.UT.f, qtl.UT.m)

qtl.UM.f <- scanone(females.UM, pheno.col = phe, method="hk")
qtl.UM.m <- scanone(males.UM, pheno.col = phe, method="hk")
plot(qtl.UM.f, qtl.UM.m)


plot(qtl.JL.m,qtl.UT.m, qtl.UM.m, main = paste0("QTL ",phe," (males)"))
plot(qtl.JL.f,qtl.UT.f, qtl.UM.f, main = paste0("QTL ",phe," (females)"))

lod.overview <- cbind(JL.M = qtl.JL.m$lod, UT.M = qtl.UT.m$lod, UM.M = qtl.UM.m$lod,
                      JL.F = qtl.JL.f$lod, UT.F = qtl.UT.f$lod, UM.F = qtl.UM.f$lod)
cor(lod.overview)

addcovar <- cbind(as.numeric(pull.pheno(mcross)[, "Site"]), as.numeric(pull.pheno(mcross)[, "Cohort.Year"]))

res <- scanone(mcross, pheno.col=phe,addcovar=addcovar, method="hk")


heatmap(t(lod.overview), scale = "none", Colv=NA)
heatmap(t(lod.overview), scale = "none", Colv=NA, breaks = c(0,3,1000), col=c("white", "black"))
