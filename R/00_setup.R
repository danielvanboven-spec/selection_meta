# ======================================================================
# 00_setup.R - PACKAGES, PATHS, GLOBAL OPTIONS
#
# Sourced first by run_all.R. Loads every package used anywhere in the
# pipeline, creates the output directories, and fixes the RNG seed so
# that all bootstraps and MCMC runs are reproducible.
#
# Produces: PATHS, SEED, N_BOOT, R_WITHIN, DOUBLE_GAMMA_II
# ======================================================================


# -- 1. Packages -------------------------------------------------------

pkgs_core <- c(
  "readxl", "readr", "dplyr", "tidyr", "tibble", "Matrix"
)
pkgs_model <- c(
  "metafor", "clubSandwich", "MCMCglmm", "brms", "posterior", "lme4"
)
pkgs_phylo <- c(
  "rotl", "ape", "ggtree", "rphylopic", "ggimage"
)
pkgs_plot <- c(
  "ggplot2", "patchwork", "cowplot", "ggbeeswarm", "orchaRd", "magick"
)

invisible(lapply(c(pkgs_core, pkgs_model), function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Missing required package: ", p, call. = FALSE)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

# Phylogeny and plotting packages are only needed by 03_phylogeny.R and
# the figures/ scripts. Load them lazily so the modelling pipeline runs
# on a machine without ggtree / rphylopic installed.
load_optional <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    warning("Optional packages not installed: ", paste(missing, collapse = ", "),
            "\n  These are needed only for the phylogeny and figure scripts.",
            call. = FALSE)
    return(FALSE)
  }
  invisible(lapply(pkgs, function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))))
  TRUE
}


# -- 2. Paths ----------------------------------------------------------

PATHS <- list(
  data       = "data",
  masterdata = "data/Masterdata.xlsx",
  raw_studies = "data/raw_studies",
  out        = "output",
  figures    = "output/figures",
  tables     = "output/tables",
  models     = "output/models"
)

for (p in PATHS[c("out", "figures", "tables", "models")])
  dir.create(p, showWarnings = FALSE, recursive = TRUE)


# -- 3. Analysis constants ---------------------------------------------
# Every tunable number the manuscript quotes lives here, so that the
# reported value and the value used in the code cannot drift apart.

SEED   <- 42          # global RNG seed
N_BOOT <- 1000        # bootstrap replicates. Methods: "B = 1,000 replicates"
                      # for the folded-normal CI and the parametric gamma
                      # bootstrap in the Krzanowski analysis.

R_WITHIN <- 0.5       # within-pair sampling correlation between gamma_ij and
                      # z(Rp) in the bivariate model. Methods: "I set
                      # r_w = 0.5 for the primary analysis".
R_WITHIN_SENS <- c(0, 0.3, 0.5, 0.7, 0.9)   # sensitivity grid

EXCLUDED_STUDY_ID <- 8  # study dropped for implausibly large gradients


