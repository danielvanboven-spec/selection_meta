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
└── output/
   └─── figures/  tables/  models/  sessionInfo.txt
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

## Dependencies

Core: `readxl`, `readr`, `dplyr`, `tidyr`, `tibble`, `Matrix`.
Modelling: `metafor`, `clubSandwich`, `MCMCglmm`, `brms`, `posterior`, `lme4`.
Phylogeny: `rotl`, `ape`, `ggtree`, `rphylopic`, `ggimage`.
Figures: `ggplot2`, `patchwork`, `cowplot`, `ggbeeswarm`, `orchaRd`, `magick`.

`brms` needs a working Stan toolchain. `ggtree` comes from Bioconductor.
`orchaRd` comes from GitHub (`daniel1noble/orchaRd`).

Analyses were run under R 4.4.1.

---

=======
# Correlational selection is widespread across plants and animals and shapes phenotypic trait covariation

Data and R code for a systematic review and meta-analysis of non-linear selection gradients, estimated within the Lande-Arnold framework.

**Authors:** Daniël van Boven, Anna Penna and Erik I. Svensson

**Contact:** erik.svensson@biol.lu.se

**Status:** manuscript in preparation. The numbers produced by this code may still change until the manuscript is submitted.

## What the analysis does

The project has three aims:

1. Estimate the magnitude of correlational selection gradients (γij) and test what predicts it, including the type of selection (natural or sexual).
2. Test whether the sign of γij matches the sign of the phenotypic correlation (Rp) between the same two traits.
3. Measure the geometric alignment between the selection surface (γ) and phenotypic trait space (P), using Krzanowski subspace comparison.

## Repository structure

| Folder | Contents |
|---|---|
| `data/` | `Masterdata.xlsx`, the full dataset |
| `R/` | All analysis and figure scripts |
| `output/figures/` | Figures as produced by the scripts |

Tables and fitted model objects are written to `output/tables/` and `output/models/` when the scripts run. Those folders are not tracked in this repository.

## The data

`data/Masterdata.xlsx` has one row per selection gradient. A row is either a quadratic gradient (γii, one trait) or a correlational gradient (γij, two traits), taken from the published study. Each row also records:

- taxonomy of the study species (kingdom, class, order, Latin and English name)
- study design: sample size, number of model parameters, setting (field, laboratory or artificial), and the grouping the original authors used (for example sex, population or year)
- the fitness measure, its type (reproductive, survival or performance), and the type of selection (natural or sexual)
- the traits, their functional type, and their means and standard deviations
- the gradient and its standard error
- elements of the phenotypic covariance matrix (P), the phenotypic correlation (Rp) and, where reported, the genetic covariance matrix (G)
- a study ID and the source

Two exclusions are applied in `R/01_load_clean.R`, not in the spreadsheet: one study with implausibly large gradients is dropped, and so are rows with a gradient outside the range -2 to 2.

## Running the analysis

The scripts use paths relative to the repository root (for example `data/Masterdata.xlsx`), so set the working directory there first.

1. Open R with the repository root as working directory.
2. Source `R/00_setup.R` and `R/00_functions.R`. These load the packages, set the paths and define the shared functions, colour palettes and constants.
3. Run the numbered scripts in `R/` in numerical order. Each script lists at the top which objects it needs and which earlier script creates them.
4. Run the figure scripts (`fig...R`) last. They read the objects created by the numbered scripts.

Notes:

- The phylogeny is built from the Open Tree of Life through the `rotl` package, so that step needs an internet connection.
- The Krzanowski angles use a parametric bootstrap and the Bayesian models use MCMC. Both take several minutes.

## Software

The analysis was run in R 4.4. Main packages:

- **Meta-analysis:** `metafor`, `MCMCglmm`, `brms`, `clubSandwich`, `orchaRd`
- **Phylogeny:** `rotl`, `ape`, `ggtree`, `rphylopic`
- **Data handling:** `dplyr`, `tibble`, `readxl`
- **Figures:** `ggplot2`, `ggbeeswarm`, `cowplot`, `patchwork`, `magick`

## Reuse

No licence has been set yet. Please contact the authors before reusing the data or code.
>>>>>>> 96ee146aaed6610a8f8bbc929147f32538e70455
