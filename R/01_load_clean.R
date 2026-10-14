# ======================================================================
# 01_load_clean.R - LOAD AND CLEAN MASTERDATA
#
# Produces: df - cleaned, harmonised, with derived variables
#
# Manuscript sections served:
#   Methods > Data extraction
#   Methods > Effect size calculation (the |gamma| <= 2 exclusion)
# ======================================================================


# -- 1. Load -----------------------------------------------------------

Masterdata <- readxl::read_excel("data/Masterdata.xlsx", sheet = 1)
df <- data.frame(Masterdata)

df <- df[, !grepl("^X\\.*[0-9]*$", names(df)), drop = FALSE]

n_raw <- nrow(df)


# -- 2. Type coercion --------------------------------------------------

numeric_cols <- c(
  "Gamma.ii", "Gamma.ij", "Gamma.ii.se", "Gamma.ij.se",
  "Pii", "Pii.SE", "Pij", "Pij.SE",
  "Gii", "Gij", "Rp", "Significance.Rp", "SE.Rp",
  "Mean.trait.1", "Mean.trait.2", "SD.trait.1", "SD.trait.2",
  "Sample.size", "No..of.parameters"
)
numeric_cols <- intersect(numeric_cols, names(df))
df[numeric_cols] <- lapply(df[numeric_cols],
                           function(x) suppressWarnings(as.numeric(as.character(x))))


# -- 3. Taxonomic harmonisation ----------------------------------------
# Only spelling is corrected here, so df carries the species that were
# actually studied and the species count is the true one. Species the
# Open Tree of Life cannot place get a stand-in in 03_phylogeny.R, for
# the position of their tip only.

df$Latin.name <- trimws(df$Latin.name)

df$Latin.name <- dplyr::recode(df$Latin.name,
                               "Ramphomyia sulcata" = "Rhamphomyia sulcata",
                               "Cypripedium acuale" = "Cypripedium acaule"
)


# -- 4. Exclusions -----------------------------------------------------
# Results: "One study was subsequently excluded because the calculated
# effect sizes were implausibly large (gamma_ij >= 2)".

df <- df[!is.na(df$ID) & df$ID != EXCLUDED_STUDY_ID, ]
n_after_study_drop <- nrow(df)

# Row-level filter: keep a row only if its gradients lie within -2 to 2.
# A missing gradient passes, so quadratic rows are not dropped for
# having no gamma_ij and the other way round.
gamma_limit <- 2

keep_ii <- is.na(df$Gamma.ii) | abs(df$Gamma.ii) <= gamma_limit
keep_ij <- is.na(df$Gamma.ij) | abs(df$Gamma.ij) <= gamma_limit

df_gamma_dropped <- df[!(keep_ii & keep_ij), ]
df <- df[keep_ii & keep_ij, ]
n_after_gamma_filter <- nrow(df)

# -- 5. Derived variables ----------------------------------------------

# "N/A" strings to real NA in the trait-2 column, so that a row with no
# second trait is unambiguously a quadratic (gamma_ii) row
df$Trait.2[df$Trait.2 == "N/A"] <- NA

# Sign of the correlational gradient
df$sign_ij <- ifelse(df$Gamma.ij > 0, 1, 0)

# Canonical trait-type pair (A x B == B x A)
df <- df %>%
  dplyr::mutate(
    trait1 = pmin(Type.of.trait, Type.of.trait2),
    trait2 = pmax(Type.of.trait, Type.of.trait2),
    pair   = paste(trait1, trait2, sep = " - ")
  )

# Type match: are the two traits in a pair of the same functional type?
# Results: "650/1010 same, 64.4%; 360/1010 different, 35.6%"
df$Type.match <- ifelse(
  is.na(df$Type.of.trait2) | df$Type.of.trait2 == "N/A", NA_character_,
  ifelse(df$Type.of.trait == df$Type.of.trait2, "Same", "Different")
)

# Within-study grouping level: replicate, sex, population, treatment or
# other stratum used by the original authors. Used as a random effect.
df$group_id <- ifelse(
  is.na(df$Grouping) | df$Grouping == "N/A",
  paste0("study_", df$ID, "_all"),
  paste0("study_", df$ID, "_", df$Grouping)
)


# -- 6. Report ---------------------------------------------------------

cat("\n---- 01_load_clean ----\n")
cat("Masterdata rows read        :", n_raw, "\n")
cat("After dropping study ID", EXCLUDED_STUDY_ID, "  :", n_after_study_drop, "\n")
cat("After |gamma| <=", gamma_limit, "filter    :", n_after_gamma_filter, "\n")
cat("  rows dropped by filter    :", nrow(df_gamma_dropped), "\n")
cat("  gamma_ii estimates        :", sum(!is.na(df$Gamma.ii)), "\n")
cat("  gamma_ij estimates        :", sum(!is.na(df$Gamma.ij)), "\n")
cat("  studies                   :", dplyr::n_distinct(df$ID), "\n")
cat("  species (Latin.name)      :", dplyr::n_distinct(df$Latin.name[!is.na(df$Latin.name)]), "\n")

cat("\nRows dropped by the gamma filter:\n")
print(df_gamma_dropped[, c("ID", "Grouping", "Trait.1", "Trait.2",
                           "Gamma.ii", "Gamma.ij")], row.names = FALSE)
