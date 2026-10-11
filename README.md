# Correlational selection and the structure of the fitness landscape

Analysis pipeline for a systematic review and meta-analysis of non-linear
selection gradients within the Lande and Arnold quantitative genetics
framework. The dataset spans 1,649 gradients from 70 studies, 63 species
and 13 taxonomic classes.

Written in R. Manuscript prepared in LaTeX, reported following
PRISMA-EcoEvo.

---

## Quick start

```r
# from the project root
source("run_all.R")
```

`run_all.R` has three switches at the top. Set `RUN_SLOW <- FALSE` to skip
the bootstraps, the MCMC chains and the raw-data refits; the pipeline then
finishes in under a minute and still produces every number that does not
depend on them.

Required: `data/Masterdata.xlsx`.
Optional: `data/raw_studies/` (script 04) and `data/photos/` (Figure 1).

---

## What changed in this restructure

The previous 25 scripts have been reorganised into 15 analysis scripts and
4 figure scripts. **Only code that produces something appearing in the
manuscript or the supplement has been kept.** Nothing has been added to the
analysis; the changes are pruning, ordering, bug fixes, and the addition of
a reconciliation step.

Read the two sections below before running anything for the submission.

---

## Layout

```
selection_meta/
├── run_all.R                  runs everything in order
├── README.md
├── R/
│   ├── 00_setup.R             packages, paths, constants (SEED, N_BOOT, R_WITHIN)
│   ├── 00_functions.R         all shared helpers
│   ├── 01_load_clean.R        df
│   ├── 02_sample_descriptives.R
│   ├── 03_phylogeny.R         tree, phylo_cor_matrix
│   ├── 04_lande_arnold_raw.R  raw-data refits, sampling covariance M
│   ├── 05_build_combined.R    combined_df, V_combined
│   ├── 06_meta_gamma_ij.R     AIM 1
│   ├── 07_gamma_ii.R          quadratic benchmark
│   ├── 08_bias_correction.R   folded-normal E[|gamma|]
│   ├── 09_publication_bias.R  Egger regressions
│   ├── 10_rp_metaregression.R AIM 2, plus sign concordance
│   ├── 11_rp_bivariate.R      AIM 2, rho_B by REML
│   ├── 12_rp_bivariate_brms.R AIM 2, rho_B by Stan
│   ├── 13_krzanowski_angles.R AIM 3, all_angles
│   ├── 14_alignment_models.R  AIM 3, angle and similarity models
│   └── 15_reconciliation.R    checks every reported number
├── figures/
│   ├── fig01_phylogeny.R
│   ├── fig02_orchard_quadrant.R
│   ├── fig03_krzanowski.R
│   └── figS4_funnel_plots.R
├── data/
│   ├── Masterdata.xlsx
│   ├── raw_studies/           the ten raw-data study files
│   └── photos/                Genus_species.jpg for Figure 1
├── output/
│   ├── figures/  tables/  models/  sessionInfo.txt
└── archive/
    └── the removed scripts, with a note on why each went
```

---

## Where each reported number comes from

### Results, sample and study characteristics

| Number | Script |
|---|---|
| 1,649 gradients; 639 gamma_ii, 1,010 gamma_ij; 70 studies | `01_load_clean.R` |
| 63 species, 13 classes, 38 orders, 2 kingdoms | `02_sample_descriptives.R` |
| 41 animal and 29 plant studies | `02_sample_descriptives.R` |
| 213 gamma-matrices, 106 full rank, 75 with P, 24 with k = 2 | `02_sample_descriptives.R`, `13_krzanowski_angles.R` |
| setting, fitness, selection, trait type, type match percentages | `02_sample_descriptives.R` |

### Results, strength and generality (aim 1)

| Number | Script |
|---|---|
| E\[\|gamma_ij\|\] = 0.087 \[0.036, 0.143\] natural | `08_bias_correction.R` |
| E\[\|gamma_ij\|\] = 0.234 \[0.063, 0.332\] sexual | `08_bias_correction.R` |
| raw \|gamma_ij\| = 0.115 \[0.088, 0.143\] and 0.217 \[0.122, 0.312\] | `06_meta_gamma_ij.R` |
| backward elimination retains type of selection only | `06_meta_gamma_ij.R` |
| QE(1008) = 8917.8; 42 / 37 / 20 per cent variance | `06_meta_gamma_ij.R` |

### Results, bivariate model (aim 2)

| Number | Script |
|---|---|
| 535 trait pairs from 38 studies | `10_rp_metaregression.R` |
| signed gamma_ij = 0.013, z(Rp) = 0.142 | `11_rp_bivariate.R` |
| REML rho_B = 0.78, weakly identified | `11_rp_bivariate.R` |
| Bayesian rho_B = 0.63 \[-0.07, 0.94\], P(rho > 0) = 0.97 | `12_rp_bivariate_brms.R` |
| slope 0.113, P = 0.017; with log n, 0.131, P = 0.005 | `10_rp_metaregression.R` |
| MCMCglmm slope 0.115, pMCMC = 0.01 | `10_rp_metaregression.R` |
| 55 per cent concordant, chi2 = 5.52, P = 0.019 | `10_rp_metaregression.R` |

### Results, angles (aim 3)

| Number | Script |
|---|---|
| 75 comparisons from 26 studies | `13_krzanowski_angles.R` |
| intercept 48.5 deg, slope -20.2 deg per unit eigenvalue | `14_alignment_models.R` |
| sigma2 203.3 and 281.2, QE(73) = 5185.7 | `14_alignment_models.R` |
| CR2 t(20.1) = 10.01 and t(3.02) = -4.67, P = 0.018 | `14_alignment_models.R` |
| similarity -0.81 logit, 0.31 back-transformed, slope 0.83 | `14_alignment_models.R` |
| mean gain 45.7 deg at k >= 2 | `13_krzanowski_angles.R`, `fig03` |

### Supplement

| Item | Script |
|---|---|
| Table S1, benchmark papers | not code; compiled by hand |
| Figure S2, PRISMA diagram | not code; drawn in LaTeX |
| Table S3, predictor variables | not code; `02_sample_descriptives.R` gives the counts |
| Figure S4, funnel plots | `09_publication_bias.R` and `figS4_funnel_plots.R` |

---

## Known discrepancies

Nine items. The first four are gaps between the manuscript and the current
`Masterdata.xlsx`, and one side of each is stale. The remaining five are
decisions in the code that the Methods describe differently. Work out which
side is wrong before editing either.

**1. Total effect sizes and study count.** The Abstract and Results say
"1,701 effect sizes from 71 studies". The current file holds 1,681 rows
carrying an ID, across 73 IDs. Dropping study 8 leaves 1,660 rows and 72
studies, and the `|gamma| <= 2` filter then empties IDs 12 and 39, giving
1,649 rows from 70 studies. That final figure matches the manuscript
exactly; the 1,701 and the 71 do not correspond to anything in the file as
it now stands. Most likely the file has been edited since the sentence was
written.

**2. Species count.** `df` carries 65 distinct binomials after the
taxonomic recoding. The manuscript says 63, and Figure 1 shows 63 tips. The
gap opens in `03_phylogeny.R`: two species survive `tnrs_match_names()` but
fail `is_in_tree()` and are dropped from the induced subtree. So 63 is the
count of species in the phylogeny, not in the dataset. Either report 65 in
the Results and 63 in the caption, or say why two were excluded. The script
prints which two.

**3. Gamma-matrix counts.** The Methods say 213 gamma-matrices, of which
106 were full rank. Recomputing on the current file with the grouping key
in `02_sample_descriptives.R` gives 176 and 119. The key is
study x grouping x species x fitness measure; dropping species from the key
gives 170. Neither reaches 213. This needs re-deriving from whatever key
produced the original numbers.

**4. Trait pairs in the Rp analyses.** The Results say 535 pairs from 38
studies. The filter in `10_rp_metaregression.R` returns 516 pairs from 38
studies on the current file. The study count matches; the pair count is 19
short.

**5. Doubling gamma_ii for the Krzanowski comparison.** The Methods state
that gamma_ii values were doubled before being placed on the diagonal, so
that quadratic and correlational terms share a scale. The original
`08_meta_curvature.R` did not double them. Controlled by `DOUBLE_GAMMA_II`
in `00_setup.R`, currently `FALSE`, which reproduces the numbers in the
manuscript. Note that the four raw-data studies already have their
quadratic coefficients doubled at the point of estimation, so switching
this on would double those twice; if you do switch it on, exclude the
recalculated studies from the second doubling.

**6. Which model draws Figure 2b.** The panel's fitted line and its slope
annotation come from `rma_rp_z_con`, which includes log sample size and
uses the conservative sampling variance. The Results attribute
slope = 0.113, P = 0.017 to the gamma-only model. The two slopes are close,
so the label is not wrong, but the caption and the Results should name the
same model.

**7. The conservative sampling variance.** `Vz_con` was `Vz * sqrt(1.5)`. A
design effect DE = 1 + (M - 1) x ICC inflates the variance by DE, so
`Vz * 1.5` is what the stated reasoning implies; `sqrt(1.5)` inflates the
standard error by that factor instead. Controlled by
`DESIGN_EFFECT_ON_VARIANCE` in `10_rp_metaregression.R`, currently `FALSE`
to reproduce the manuscript.

**8. Four raw-data studies or ten.** The Methods say gradients were
recalculated for four studies. `04_lande_arnold_raw.R` refits ten:
4, 6, 21, 23, 25, 29, 31, 33, 40, 43. Both can be true if only four lacked
published gradients and the other six were refitted solely to obtain the
off-diagonal sampling covariances for V. If so, the Methods sentence needs
that clause. Check the `Report/Calculated` column for those ten IDs.

**9. Bootstrap replicates.** The Methods say B = 1,000. The bias-correction
script ran B = 100. It is now driven by `N_BOOT` in `00_setup.R`, so the
bootstrap CIs will shift slightly on the first rerun at 1,000. Expect the
intervals to tighten a little; the point estimates are unaffected.

---

## Dependencies

Core: `readxl`, `readr`, `dplyr`, `tidyr`, `tibble`, `Matrix`.
Modelling: `metafor`, `clubSandwich`, `MCMCglmm`, `brms`, `posterior`, `lme4`.
Phylogeny: `rotl`, `ape`, `ggtree`, `rphylopic`, `ggimage`.
Figures: `ggplot2`, `patchwork`, `cowplot`, `ggbeeswarm`, `orchaRd`, `magick`.

`brms` needs a working Stan toolchain. `ggtree` comes from Bioconductor.
`orchaRd` comes from GitHub (`daniel1noble/orchaRd`).

Analyses were run under R 4.4.1.

---

