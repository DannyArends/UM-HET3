## UM-HET3 ✨
Data analysis UM-HET3 - 4-way cross created from (BALB/cByJ x C57BL/6J) and (C3H/HeJ x DBA/2J)

Creating a physical and genetic map from Monsterplex Capture DNA-Sequencing data and Sequenom MassARRAY data for the UM-HET3 populations

- Three physical locations (Sites/Centers): JL, UM, and UT
- Cohort Structure 2004 - 2012
- Treatments (21 levels, some duplicated/too small)
- Body weight at 42 days, and 6, 12, 18, and 24 months 
- Individual lifespan data in days

Progressive QTL mapping code can be found in this repository

Genotype coding convention in R/1tl v1: BALB/cByJ = AA, C57BL/6J = BB, C3H/HeJ = CC, DBA/2J = DD

Phase-known alleles, with uncertain/not imputed counterpart:
AC or AD = A?
BC or BD = B?
AC or BC = ?C
AD or BD = ?D

### Structure 📁

Most scripts in this repository start by loading in the um-het3-rqtl.csvr, which is the data coded for R/qtl v1 as explained in the read.cross help file. Data can be loaded in with read.cross parameters set to genotypes=NULL and na.strings=c("-", "NA")

The folder [PreCross/](./PreCross/) contains all the code used to call SNPs on the monsterplex BAM files, and the conversion & harmonization of Monsterplex and Sequenom data into the cross object. It also contains the code to convert observed SNPs to the Founder strain haplotypes.

The folder [ProgessiveMapping/](./ProgessiveMapping/) contains all the code used to perform actuarial QTL scans on lifespan for 4-way, paternal and maternal maps.

The folder [Soma/](./Soma/) contains all the code used to perform the Correlated Trait Locus (CTL) mapping analysis between bodyweight (at 5 timepoints) and progressive lifespan mapping.

The folder [Inversions/](./Inversions/) contains the inversion finder code, which finds inversion by looking at recombination frequency using an 8-marker sliding window

The folder [Old/](./Old/) contains the old & deprecated code not used anymore, or which was used to do one off analysis.

### Contributing 🙌

Want to contribute? Great! Contribute to this repo by starring ⭐ or forking 🍴, and feel 
free to start an issue first to discuss idea's before sending a pull request. You're also 
welcome to post comments on commits.

Or be a maintainer, and adopt (the documentation of) a function.

### License ⚖️

Written by Danny Arends and is released under the GNU GENERAL PUBLIC LICENSE Version 3 (GPLv3). 
See [LICENSE.txt](./LICENSE.txt).
