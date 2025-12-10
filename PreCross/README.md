## PreCross ✨

This folder contains all the code used to call SNPs on the monsterplex BAM files, and the conversion & harmonization of Monsterplex and Sequenom data into the cross object. It also contains the code to convert observed SNPs to the Founder strain haplotypes. This code is used to produce the um-het3-rqtl.csvr from Monsterplex .BAM files and the Sequenom data, with phenotypes and covariates coming from
[genenetwork.org](https://genenetwork.org).

Since we annot easily distribute BAM files in a github repository, please contact Rob W Williams <labwilliams@gmail.com> or Danny Arends <Danny.Arends@gmail.com> to get a copy of the input data. Code in this folder is not easily generalizable, since it depends on mouse genome project files as well as certain reference genome files to be available and indexed at specific location.

### Contributing 🙌

Want to contribute? Great! Contribute to this repo by starring ⭐ or forking 🍴, and feel 
free to start an issue first to discuss idea's before sending a pull request. You're also 
welcome to post comments on commits.

Or be a maintainer, and adopt (the documentation of) a function.

### License ⚖️

Written by Danny Arends and is released under the GNU GENERAL PUBLIC LICENSE Version 3 (GPLv3). 
See [LICENSE.txt](./LICENSE.txt).

### Cite 📝

Arends, D., Ashbrook, D. G., Roy, S., Lu, L., Sloan, Z., Centeno, A. G., Lamour, K. H., de Magalhães, J. P., Prins, P., Broman, K. W., Sen, S., Mitchell, S. J., MacArthur, M. R., Akin, Ö. A., Auwerx, J., Bajwa, A., Diaz, V., Harrison, D. E., Strong, R., Nelson, J. F., Mozhui, K., Williams, E. G., Miller, R. A., & Williams, R. W. (2025). **Genetic Modulation of Lifespan: Dynamic Effects, Sex Differences, and Body Weight Trade-offs** *bioRxiv* [https://doi.org/10.1101/2025.04.27.649857](https://doi.org/10.1101/2025.04.27.649857)

