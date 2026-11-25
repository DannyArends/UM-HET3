# Errata (draft): Arends et al. 2026, *Nature* 654:437

Draft errata for *Dynamics of genetic and somatic trade-offs in ageing and mortality* ([doi.org/10.1038/s41586-026-10407-9](https://doi.org/10.1038/s41586-026-10407-9)).

**Status: draft of 25 November 2026, prepared by Claude (Opus 5.5) for Rob Williams. Not yet reviewed by the other authors; nothing here has been sent to *Nature*.**

- **Sources:** the published PDF, Supplementary Tables 7, 9 and 12, `DataSet/um-het3-rqtl.csvr` and `DataSet/regions/` in this repository, and reruns of the published analyses.
- **Confidence:** **Sure** (two sources in the paper or its data disagree and one is right), **Likely** (one reading clearly fits), **Ask** (an author should decide).
- **Reruns behind section A:** [Cauchy review](https://github.com/robwwilliams/bxd_mortality_rates_2026/blob/main/notes/Cauchy_review_2026-11-21.md) and [code](https://github.com/robwwilliams/bxd_mortality_rates_2026/blob/main/UMHET3/cauchy_review.R); [permutation results](https://github.com/robwwilliams/bxd_mortality_rates_2026/blob/main/output/UMHET3_permutations.txt) and [code](https://github.com/robwwilliams/bxd_mortality_rates_2026/blob/main/UMHET3/umhet3_permutations.R); [Chr 1 genotype calls](https://github.com/robwwilliams/bxd_mortality_rates_2026/blob/main/notes/Message_Chr1_genotype_error_2026-10-26.md). The permutations are labeled exploratory.

## A. Statistics and data: corrections of substance

| # | Where | Now | Should be, and why | Confidence |
|---|---|---|---|---|
| A1 | Methods, "Significance testing"; Table 1 (Cauchy column); `ActuarialMapping/Cauchy.R` | LODs enter the Cauchy test as *P* = 10^−LOD | These LODs are on 3 d.f., so *P* = Pr(χ²₃ > 2 ln10 · LOD). Cauchy −log*P* falls by 0.45–0.72 at every locus at or above 1.5. The number at ≥ 1.5: combined 25 → 22, male 25 → 22, female 9 → 5. **Vita3a (best 1.39, male) and Vita11b (best 1.43, female) fall below 1.5 in every column**; ten locus-by-sex results fall in all (see the table below). | Sure |
| A2 | Methods; Table 1 footnote e; Table 2 footnote b; Fig. 3b legend | "LOD 3.65, 3.95 and 4.65" for genome-wide *P* 0.10, 0.05, 0.01 | These are −log10 *P* values over 442 tests, applied as LODs. The scheme in `ActuarialMapping/stratifiedPerms.R` (1,000 shuffles within sex × site, run on the published cross) gives **3.73–3.74, 4.08–4.12 and 4.86–4.97** for one T-age. The Bonferroni version on 3 d.f. would be 4.21, 4.53, 5.26, which is stricter than the permutations. Suggested: the permutation values. | Sure (numbers); Ask (which to adopt) |
| A3 | Table 1 (peak LODs) and the text "29 Vita loci" | All 29 presented as significant | At the peak LOD, T-age and first-listed sex in Table 1, genome-wide *P* from the permutations is above 0.05 for **Vita5a (0.072) and Vita11b (0.068)**. Vita3a is at 0.048 with the Table 1 LOD of 4.26, but its largest LOD in Supplementary Table 9 is 3.83 (P = 0.11; see question 4), which would not pass. So Vita3a, 5a and 11b are the loci to settle. **The "any T-age" test is stricter:** at 5% it keeps 12 of 29 loci and at 10%, 17 (best LOD in the locus interval, any map; where intervals overlap, a locus can borrow its neighbor’s peak, so these are upper limits). | Sure |
| A4 | Table 1 Cauchy column, male | Vita14a and Vita18a are 0.00 | A LOD rounded to 0.00 sets the combined *P* to 1. With that T-age left out: 3.08 and 2.34 (2.48 and 1.81 on 3 d.f.). | Sure |
| A5 | Chr 1 genotypes: Vita1a, 1b, 1c?, 1d; Supplementary Tables 2, 5, 7 | Paternal H : D is about 2:1 at Vita1a and Vita1b, 1.4:1 at Vita1d; almost no BD at Vita1b | One of two marker sets (set B) loses paternal-D heterozygotes at SNPs 1_3010274, 1_31062737, 1_132295971, 1_167148678 and 17_83900092, so BD is read as BH. UM and UT are affected, JL is not. A rescan of Vita1a–Vita1d without the faulty calls is needed; direction of change unknown. Detail: [Chr 1 genotype calls](https://github.com/robwwilliams/bxd_mortality_rates_2026/blob/main/notes/Message_Chr1_genotype_error_2026-10-26.md). | Sure (miscalls); Ask (effect on LODs) |

**A1 and A3 together, the loci that lose support:**

| Locus | Best published Cauchy → on 3 d.f. | Permutation *P* at Table 1's peak (one T-age) | Any T-age |
|---|---|---|---|
| **Vita11b** | 1.96 → 1.43 (F); fails in all columns | 0.068 | 0.58 |
| **Vita3a** | 1.91 → 1.39 (M); fails in all columns | 0.048 (0.11 with the LOD of 3.83) | 0.54 |
| **Vita5a** | 2.10 → 1.56 (M) | 0.072 | 0.68 |
| Vita6a | 2.58 → 2.02 (M) | 0.040 | 0.62 |
| Vita11a | 2.21 → 1.67 (M) | 0.038 | 0.57 |
| Vita10a | 2.48 → 1.92 (C) | 0.028 | 0.36 |
| Vita4b | 2.38 → 1.84 (M) | 0.025 | 0.44 |

*Results by sex that drop below 1.5 on the corrected Cauchy:* combined Vita11a, 11c, 15b; male Vita2c, 3a, 9a; female Vita1a, 9b, 11b, 12a. Most of these loci still pass in another column.

## B. Text, legends and tables that state something wrong

| # | Where | Now | Should be | Evidence | Confidence |
|---|---|---|---|---|---|
| B1 | Fig. 1h,i legend | "maternal block with mortality rates of the **C and D** maternal haplotypes" | C and **B** | UM-HET3 mothers give C or B; Extended Data Fig. 4 says "C and B" | Sure |
| B2 | Table 4, Soma14b row | Vita overlap "Vita14**a**" | Vita14**b** | Same marker, 14_101437466 (rs49624430); Supplementary Table 12 says Soma14b = Vita14b. *An off-by-one letter* | Sure |
| B3 | Table 2, Vita6b, males | Partner "15**a**" | "15**b**" | Male LOD Vita6b × 15b 4.21, × 15a 2.85 (Supp. 12.9); Vita15b's own row lists 6b. *Off by one* | Sure |
| B4 | Table 2, Vita11a, females | "–" | "1a, 14a" | LOD 11a × 1a 4.68, × 14a 4.26; Vita1a's row and the text list 11a | Sure |
| B5 | Table 2, Vita11b, females | "1a, 14a" | "6a, 14a" | LOD 11b × 1a 1.40, × 6a 4.03, × 14a 4.77 | Sure |
| B6 | Table 2, Vita11b, males | "–" | "17a" | LOD 3.92; Vita17a's row lists 11b | Sure |
| B7 | Table 2, males: omitted pairs | Not listed | 9b × 2c (4.02) and 14a × 2b (4.24) | Both reach the footnote's LOD 3.8 | Sure |
| B8 | Table 2, Vita5a, females | "3a, 17" | "3a, 17**a**" | Vita17a's row lists 5a | Sure |
| B9 | Tables 3 and 4, Soma6a | Chr **5** | Chr **6** | Its marker is 6_8006720 (and the name is Soma6a) | Sure |
| B10 | Fig. 5f panels and legend | Panels "T740", LOD 5.6 (female) and 0.9 (male); legend "weak in the T42 female survivorship (e) but strong at T740 (f)" | Either the panels are T42, or the LODs are wrong | Vita3a × Vita5a: female 5.64 at T42 and 2.72 at T740; male 0.86 and 3.14. So the numbers shown are the T42 ones and the interaction is strong at T42, the opposite of the legend | Sure |
| B11 | Fig. 4i legend | "not in the older **T665** survivorship" | T**865** (as in panel n and the next sentence) | Same legend, panels m,n | Sure |
| B12 | Fig. 4 legend, d–h | "42, 183, **385**, **558** and 730 days" | 365 and 548 | Methods and main text; panel g says 548 | Sure |
| B13 | Fig. 4 legend | "**LOD** scores are significant at 2.75" | –log*P* | Table 3 footnote; a CTL test, not LOD | Likely |
| B14 | Table 1, Vita9c (T980, combined) | *n* = 1,990, mean 965 days | 1,812 and 1,074 | Mean of mice alive at T980 cannot be below 980; *n* and mean from Supplementary Table 7 and `um-het3-rqtl.csvr` | Sure |
| B15 | Table 1, Vita9b (T1025, combined) | *n* = 1,847, mean 1,112 | 1,246 and 1,106 | More mice than there are alive at T980 (1,812) | Sure |
| B16 | Table 1, Vita15a (T905, combined) | *n* = 3,549, mean 1,023 | 2,740 and 1,030 | Same source | Sure |
| B17 | Table 1, Vita14b (T80, combined) | *n* = 3,377, mean 1,006 | 6,434 and 845 at T80 (3,377 and 1,002 are T860, the T-age of the highest combined LOD per footnote b) | Same source | Sure |
| B18 | Table 1, Vita9a (T95, female) | mean 910 (*n* = 2,754) | 888 (*n* = 3,035); 2,754 is the female *n* at T665 | Same source | Likely |
| B19 | Table 1, others | Vita2a *n* 1,548 (male T140 is 3,396); Vita1c *n* 3,255 (3,379); Vita3a 368 (383); Vita5a 350 (333); VitaXa 3,332 and 806 (3,362 and 813; 3,332 is Vita11c's T42 entry); Vita13a 6,203 (6,402); Vita18a mean 858 (845); Vita17a mean 925 (918); Vita15b mean 813 (807) | Recheck each | *n* may be restricted to genotyped mice, so some may be right | Ask |
| B20 | Table 1 and Table 2, Vita6a | Marker rs38831890 at 6_93680853, peak 107.4, interval 97.3–122.3 Mb | Recheck: the marker lies outside the interval; `DataSet/regions/regions_4way_topm.txt` gives 6_107382038 | Table 1; `DataSet/regions/` | Likely |
| B21 | Table 3: Soma8b, Soma4a | Peak 111.3 but marker 8_126505019; peak 11.2 but marker 4_30761996 | Recheck (footnote a does not mention imputed peaks) | Table 3 | Ask |
| B22 | Fig. 2e,f legend against the text | Legend: "haplotype **C** mortality starting at 400 days", reversal T**450**–T**810**. Text: "a **D** haplotype wave of mortality starts at 400 days", reversal T**410**–T**800** | Make the two agree | Both in the paper | Ask (Fig. 2f settles it) |
| B23 | Results, Loci summary | "36 ± 12 days (mean ± **s.e.m.**)" | s.d. | The 29 "Max peak diff." values: mean 35.1, s.d. 12.2, s.e.m. 2.3 | Sure |
| B24 | Results, Vita durations | "349 ± 229 days" | 328 ± 227 by summing the "Duration" column | Recheck | Ask |
| B25 | Table 2, "Fourteen Vita loci have strong G × S" | 14, with Vita9c marked X | 13 at the stated LOD ≥ 2.75; Vita9c is 2.73 | Table 2 footnote c | Likely |
| B26 | Results and Methods, Vita9b | "127 protein-coding genes with missense variants" | Methods: "262 protein-coding genes … 98 with missense variants" | Both in the paper | Ask |

## C. Cross-references that point to the wrong place

| # | Where | Now | Should be | Confidence |
|---|---|---|---|---|
| C1 | Results, Soma loci | "Fig. 5b,c", "Fig. 5i", "Fig. 5j" | Fig. **4**b,c; 4i; 4j | Sure |
| C2 | Vita1a discussion | "(Fig. 6g,h)" | Fig. **5**g,h (the circle plots) | Sure |
| C3 | Fig. 2d legend | Vita9c "(Extended Data Fig. 2**q**)" | 2**p** (2q is Vita10a). *Off by one* | Sure |
| C4 | Results, Vita2 complex | "Extended Data Fig. 1e,f" for Vita2b and Vita2c | Extended Data Fig. **2f,g** (Extended Data Fig. 1 has panels a–c; 2e,f are Vita2a and Vita2b, so this is *off by one* too) | Likely |
| C5 | Results, VitaXa | "Extended Data Figs. 1**ac** and 2ac" | Extended Data Fig. 2ac (and 3ac); there is no 1ac | Likely |
| C6 | Vita9b text (twice) | "Extended Data Fig. 2o2" | Extended Data Fig. **3**o2 (2 has single letters) | Likely |
| C7 | Extended Data Fig. 1 legend | "Other conventions as in Fig. **2d**" | Fig. **1d** | Sure |
| C8 | Table 1 footnote h | "See Extended Data Fig. **3** for variances" | Extended Data Fig. **6** | Sure |
| C9 | Dynamics of heritability | "(Fig. 2g and Extended Data Fig. **5**)" | Extended Data Fig. **6** | Likely |
| C10 | Fig. 5a legend | "for three sex combinations—combined, females and males—in Supplementary Table **10**" | Supplementary Table **12**; it holds females and males only | Sure |
| C11 | Methods, *C. elegans* | "Ortholist2 (Supplementary Table **16**)" | Supplementary Table **14** (16 is Gene Ontology) | Likely |
| C12 | Table 4 footnote e | "Vita type and sex from Table **1**" | Table **2** | Likely |
| C13 | Extended Data Fig. 5 legend | "c,d, Corresponding haplotype-specific plots" | c–f | Likely |

## D. Typographical

| # | Where | Now | Should be |
|---|---|---|---|
| D1 | Extended Data Fig. 8 legend | "are **provide** in Fig. 4i–l" | provided |
| D2 | Extended Data Fig. 5 title | "of **a** Vita14b" | of Vita14b |
| D3 | Fig. 6b,c legend | "effect plots **with for** the paternal" | for the paternal |
| D4 | Table 3 footnote b | "Max –**lopP**" | –logP |
| D5 | Table 2 footer | "**Supplentary** Information" | Supplementary |
| D6 | Vita1a discussion | "strong **constraint**" | constraints |
| D7 | Supplementary Table 12.6 title | "T**3652** survivorship" | T365 |
| D8 | Supplementary Table 12.10 tab | "**GxD**_M_365" | GxG_M_365 |
| D9 | Supplementary Table 12 notes | "Pllease", "Singificant" | Please, Significant |

## E. Questions for the authors

1. **Vita9c's start:** Fig. 2d says "late effects after T935", Table 1 starts Vita9c at T920 and Vita9b at T935. Which did the legend mean?
2. **72 survivorships:** they are in 15-day steps from T35 to T1,100 in the data and code, but the text starts at T42, which is not on that grid. Say how T42 and T35 relate.
3. **Table 2, Vita18a:** the layout may have put its female partner in the male column; one value, "1d", should be in the female column (LOD 18a × 1d is below 3.8 in males, 2.38).
4. **Peak LODs:** some in Table 1 exceed the largest LOD in Supplementary Table 9 anywhere in the locus interval (combined, male or female): Vita1a 6.84 against 6.52, Vita3a 4.26 against 3.83, Vita9b 6.66 against 5.85. If the Table 1 values are at imputed positions (footnote c), say so there.
5. **Rescan and reissue:** A1–A5 change counts and some loci. Correction, or a full amended Table 1?

## What was checked and what was not

- **Checked:** every Vita and Soma name in the running text, legends and Tables 1–4 that survives in the PDF text; every cross-reference to Extended Data and Supplementary Tables that could be found; Table 2's interactions against Supplementary Table 12 at T42; Table 1's sample sizes and means against Supplementary Table 7 (the *n* and mean in B14–B18 also against `um-het3-rqtl.csvr`); Fig. 5's labels by eye.
- **Not checked:** labels drawn as outlines in Extended Data Figs 4, 6, 7 and 9; reference list; Tables 3 and 4 numbers against Supplementary Table 10; the Supplementary Information PDF; the wording of the abstract and Methods.
