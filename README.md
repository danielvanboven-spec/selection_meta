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
