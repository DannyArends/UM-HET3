# UM-HET3
Data analysis UM-HET3 - 4way F2 between (BALB/cByJ x C57BL/6J) and (C3H/HeJ x DBA/2J)

Creating a physical and genetic map from Monsterplex Capture DNA-Sequencing data and Sequenom MassARRAY data for the UM-HET3 populations

- Three physical locations (Sites/Centers): JL, UM, and UT
- Cohort Structure 2004 - 2012
- Treatments (21 levels, some duplicated/too small)
- Body weight at 6, 12, 18, and 24 weeks and longevity data in days


Genotype coding convention: BALB/cByJ = AA, C57BL/6J = BB, C3H/HeJ = CC, DBA/2J = DD

Phase-known alleles, with uncertain/not imputed counterpart:
AC or AD = A?
BC or BD = B?
AC or BC = ?C
AD or BD = ?D

R/qtl coding as in read.cross, with genotypes=NULL and na.strings=c("-", "NA")
