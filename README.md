## Genetics of longevity

This repository contains the R code used to produce the findings in our Nature paper [doi.org/10.1038/s41586-026-10407-9](https://doi.org/10.1038/s41586-026-10407-9). 
It covers actuarial QTL mapping of *Vita* loci that shape lifespan, CTL mapping of *Soma* loci that control the body weight vs longevity trade-off, and epistatic 
interaction analysis between all loci stratified by sex and age.

Data is available in [Datasets/](./Datasets/) and on [GeneNetwork](https://genenetwork.org/show_trait?trait_id=10001&dataset=HET3-ITPPublish).

### Installation

To obtain the source code in your current working directory, run the following command your terminal:

```
git clone https://github.com/DannyArends/UM-HET3.git
```

To install R package dependencies, navigate into the UM-HET3 directory and run the following command in your terminal:

```
Rscript installer.R
```

### Structure

Most scripts in this repository start by loading in the um-het3-rqtl.csvr, which is the data coded for R/qtl-v1 as explained in the read.cross help file. 
The um-het3-rqtl.csvr file can be loaded in with the R/qtl-v1 read.cross function, with parameters set to genotypes=NULL and na.strings=c("-", "NA")

The folder [PreCross/](./PreCross/) contains all the code used to call SNPs on the monsterplex BAM files, and the conversion & harmonization of Monsterplex 
and Sequenom data into the cross object. It also contains the code to convert observed SNPs to the Founder strain haplotypes. This code is used to produce 
the um-het3-rqtl.csvr from Monsterplex .BAM files and the Sequenom data, with phenotypes and covariates coming from [genenetwork.org](https://genenetwork.org).

The folder [ActuarialMapping/](./ActuarialMapping/) contains all the code used to perform actuarial QTL scans on lifespan for 4-way, paternal and maternal 
maps used to discover *Vita* loci. This folder also contains the adjustXprobs.R code, that is needed for proper X-chromosome mapping, and furthermore 
contains bodyweight QTL mapping as well as simpleM and Cauchi Combination tests for multiple testing correction. The Lasso based non-parametric Quantile 
Regression 'engine' code is also found here.

The folder [Soma/](./Soma/) contains all the code used to perform the Correlated Trait Locus (CTL) mapping analysis between bodyweight (at 5 timepoints) 
and progressive lifespan mapping. Furthermore, code for visualizing *Soma* loci is also located in this folder.

The folder [Interactions/](./Interactions/) contains the interaction models used to investigate / map all interactions between G x G, G x Sex, G x Site, 
and G x Drug. As well as code for creating the big interaction table figure inside the paper.

The folder [Figures/](./Figures/) contains code that produces the figures / images inside the paper. Not all figure code is found here some code lives in 
their respective folders (e.g. visualizations of predictions are found in [Predictions/](./Predictions/), while some GxG images are found 
in [Interactions/](./Interactions/).

The folder [Inversions/](./Inversions/) contains the inversion finder code, which finds inversion by looking at recombination frequency using an 8-marker 
sliding window.

The folder [CandidateGenes/](./CandidateGenes/) contains all the code used to investigate and prioritize genetic features and possible candidate genes, 
code for biomaRt, VEP prediction, pathway over-representation analysis (ORA), and cross referencing genes in *Vita* / *Soma* regions with homologs in 
GenAge known to influence lifespan is found in this folder.

The folder [SingleLocus/](./SingleLocus/) contains all the code used to investigate single loci in the UM-HET3 genome (e.g. The Vita9B locus previously 
found by R. Miller), as well as code for mapping the control animals used in the Science paper by Sleiman et al. to compare/align their original work 
on a subset of animals with our mapping approach. Furthermore, there is code to investigate distal chromosome X and proximal chromosome 1, at which we 
observed some distortion in recombination frequency.

The folder [Various/](./Various/) contains all the code for computing main effect and epistatic heritability, as well code for the adaptLoess approach 
(and figures) and other various scripts that do not neatly fit into any other folder.

The folder [Predictions/](./Predictions/) contains all the code used to predict lifespan either by using random forest on covariates, as well as our 
two novel approaches ways of performing polygenic risk prediction (not in the paper, but in the biorxiv version).

The folder [Old/](./Old/) contains the old & deprecated code not used anymore, or which was used to do one off analysis.

### Dependencies

The code in this repository is made possible by, and has dependencies on, the following excellent software:

- [The R Project for Statistical Computing](https://www.r-project.org/)
- [Ensembl Variant Effect Predictor (Ensembl VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- [R/qtl](https://rqtl.org/)
- [Bcftools](https://samtools.github.io/bcftools/)

As well as many many different R-packages used in the analysis of the data.

### Genotype Coding Convention 🧬

The code for progressive QTL mapping in this repository follows specific conventions to ensure consistency and clarity. The analysis relies on a genetic 
map created from Monsterplex Capture DNA-Sequencing and Sequenom MassARRAY data. The codebase is designed to handle the complex genetic structure of the 
UM-HET3 population, a four-way cross derived from four distinct parental strains.

The data is formatted for R/qtl-v1, and as such genotypes are coded as follows:

- AA for the BALB/cByJ strain
- BB for the C57BL/6J strain
- CC for the C3H/HeJ strain
- DD for the DBA/2J strain

The code also accounts for phase-known alleles where the counterpart is uncertain or not imputed. These are represented with a question mark. 
For example, a genotype with an A allele from the BALB/cByJ parent and an unknown allele from the other is coded as A?. Similarly, an unknown allele 
combined with a known C allele is represented as ?C. This convention allows the analysis to proceed even with incomplete genotyping data, maintaining 
the integrity of the genetic map.

### Contributing
Want to contribute? Great! Contribute to this repo by starring ⭐ or forking 🍴, and feel free to start an issue first to discuss idea's before 
sending a pull request. You're also welcome to post comments on commits.

### License

Written by Danny Arends and released under the GNU GENERAL PUBLIC LICENSE Version 3 (GPLv3). 
See [LICENSE.txt](./LICENSE.txt).

### Cite

Arends, D., Ashbrook, D. G., Roy, S., Lu, L., Sloan, Z., Centeno, A. G., Lamour, K. H., de Magalhães, J. P., Prins, P., Broman, K. W., Sen, S., Mitchell, S. J., MacArthur, M. R., Akin, Ö. A., Auwerx, J., Bajwa, A., Diaz, V., Harrison, D. E., Strong, R., Nelson, J. F., Mozhui, K., Williams, E. G., Miller, R. A., & Williams, R. W. (2026). **Dynamics of genetic and somatic trade-offs in ageing and mortality** *Nature* [https://doi.org/10.1038/s41586-026-10407-9](doi.org/10.1038/s41586-026-10407-9)
