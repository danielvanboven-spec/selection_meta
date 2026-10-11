# ======================================================================
# 02_sample_descriptives.R - SAMPLE AND PREDICTOR DESCRIPTIVES
#
# Produces every count and percentage quoted in:
#   Methods > Predictor variables   (taxonomy, setting, fitness, traits)
#   Results > Sample and study characteristics
#
# Requires: df (01_load_clean.R)
# Produces: gamma_check, tax_summary, predictor_tables
#           output/tables/predictor_distributions.csv
# ======================================================================


# -- 1. Effect-size counts ---------------------------------------------

n_ii <- sum(!is.na(df$Gamma.ii))
n_ij <- sum(!is.na(df$Gamma.ij))

cat("\n---- Effect sizes ----\n")
cat("gamma_ii (quadratic)   :", n_ii, "\n")
cat("gamma_ij (correlational):", n_ij, "\n")
cat("total                  :", n_ii + n_ij, "\n")
cat("studies                :", dplyr::n_distinct(df$ID), "\n")


# -- 2. Taxonomy -------------------------------------------------------

valid <- function(x) !is.na(x) & x != "N/A"

tax_summary <- df %>%
  dplyr::summarise(
    n_species = dplyr::n_distinct(Latin.name[valid(Latin.name)]),
    n_class   = dplyr::n_distinct(Class[valid(Class)]),
    n_order   = dplyr::n_distinct(Order[valid(Order)]),
    n_kingdom = dplyr::n_distinct(Kingdom[valid(Kingdom)])
  )

cat("\n---- Taxonomic coverage ----\n")
print(as.data.frame(tax_summary), row.names = FALSE)

# Studies, species and classes per kingdom.
# Results: "41 animal studies (34 species in 11 classes) and 29 plant
# studies (29 species in 2 classes)".
tax_by_kingdom <- df %>%
  dplyr::filter(valid(Kingdom)) %>%
  dplyr::group_by(Kingdom) %>%
  dplyr::summarise(
    n_studies = dplyr::n_distinct(ID),
    n_species = dplyr::n_distinct(Latin.name[valid(Latin.name)]),
    n_class   = dplyr::n_distinct(Class[valid(Class)]),
    n_order   = dplyr::n_distinct(Order[valid(Order)]),
    .groups   = "drop"
  )

cat("\n---- By kingdom ----\n")
print(as.data.frame(tax_by_kingdom), row.names = FALSE)


# -- 3. Predictor distributions across the 1,010 gamma_ij estimates ----
# Methods: "the distributions reported below are given across the 1,010
# correlational gradients".

ij <- df[!is.na(df$Gamma.ij), ]

pct_table <- function(x, label) {
  x <- x[!is.na(x) & x != "N/A"]
  tb <- sort(table(x), decreasing = TRUE)
  data.frame(
    variable = label,
    level    = names(tb),
    n        = as.integer(tb),
    total    = length(x),
    pct      = round(100 * as.numeric(tb) / length(x), 1),
    row.names = NULL
  )
}

predictor_tables <- dplyr::bind_rows(
  pct_table(ij$Setting.field.lab.artificial, "Setting"),
  pct_table(ij$Type.of.fitness,              "Type of fitness"),
  pct_table(ij$Type.of.selection,            "Type of selection"),
  pct_table(c(ij$Type.of.trait, ij$Type.of.trait2), "Type of trait (both slots)"),
  pct_table(ij$Type.match,                   "Type match")
)

cat("\n---- Predictor distributions (gamma_ij rows) ----\n")
print(predictor_tables, row.names = FALSE)

utils::write.csv(predictor_tables,
                 file.path(PATHS$tables, "predictor_distributions.csv"),
                 row.names = FALSE)

cat("\nSample size across gamma_ij rows: median",
    stats::median(ij$Sample.size, na.rm = TRUE),
    "| range", paste(range(ij$Sample.size, na.rm = TRUE), collapse = " to "), "\n")


# -- 4. Gamma-matrix inventory -----------------------------------------
# Methods: "there were 213 gamma-matrices (defined by unique study x
# grouping x fitness-measure combinations), of which 106 were full-rank".
# A matrix is full rank when it has at least two traits, a diagonal
# entry for every trait, and every off-diagonal pair.

split_blocks <- function(cells) {
  block <- integer(length(cells)); b <- 1L; seen <- character(0)
  for (i in seq_along(cells)) {
    if (cells[i] %in% seen) { b <- b + 1L; seen <- character(0) }
    seen <- c(seen, cells[i]); block[i] <- b
  }
  block
}

df <- df %>%
  dplyr::mutate(
    cell = ifelse(is.na(Trait.2),
                  Trait.1,
                  paste(pmin(Trait.1, Trait.2), pmax(Trait.1, Trait.2), sep = " | ")),
    mkey = paste(ID, Grouping, Latin.name, Fitness.measurement, sep = "|")
  ) %>%
  dplyr::group_by(mkey) %>%
  dplyr::mutate(matrix_id = paste0(mkey, "::", split_blocks(cell))) %>%
  dplyr::ungroup()


# ---- Inventory --------------------------------------------------------

trait_count <- df %>%
  dplyr::group_by(matrix_id) %>%
  dplyr::summarise(
    traits = dplyr::n_distinct(c(Trait.1[!is.na(Trait.1)],
                                 Trait.2[!is.na(Trait.2)])),
    .groups = "drop")

diag_count <- df %>%
  dplyr::filter(is.na(Trait.2), !is.na(Gamma.ii)) %>%
  dplyr::group_by(matrix_id) %>%
  dplyr::summarise(n_diag = dplyr::n_distinct(Trait.1), .groups = "drop")

offdiag_count <- df %>%
  dplyr::filter(!is.na(Trait.2), !is.na(Gamma.ij)) %>%
  dplyr::group_by(matrix_id) %>%
  dplyr::summarise(n_offdiag = dplyr::n_distinct(cell), .groups = "drop")

# One row per matrix, so the join cannot duplicate rows
matrix_meta <- df %>%
  dplyr::group_by(matrix_id) %>%
  dplyr::summarise(
    ID                  = dplyr::first(ID),
    Grouping            = dplyr::first(Grouping),
    Latin.name          = dplyr::first(Latin.name),
    Fitness.measurement = dplyr::first(Fitness.measurement),
    P.available         = any(P.matrix.available. == "Yes", na.rm = TRUE),
    .groups = "drop")

gamma_check <- matrix_meta %>%
  dplyr::left_join(trait_count,   by = "matrix_id") %>%
  dplyr::left_join(diag_count,    by = "matrix_id") %>%
  dplyr::left_join(offdiag_count, by = "matrix_id") %>%
  dplyr::mutate(
    n_diag    = tidyr::replace_na(n_diag, 0),
    n_offdiag = tidyr::replace_na(n_offdiag, 0),
    full_rank = traits >= 2 &
      n_diag    >= traits &
      n_offdiag >= traits * (traits - 1) / 2
  )

gap_check <- dplyr::filter(gamma_check, full_rank)

cat("\n---- Gamma-matrix inventory ----\n")
cat("gamma-matrices                                     :", nrow(gamma_check), "\n")
cat("  full rank                                        :", sum(gamma_check$full_rank), "\n")
cat("  full rank AND P matrix reported                  :", sum(gap_check$P.available), "\n")
cat("  of those, >= 4 traits (eligible for k = 2)       :",
    sum(gap_check$P.available & gap_check$traits >= 4), "\n")

cat("\nP-matrix availability among full-rank matrices:\n")
print(as.data.frame(dplyr::count(gap_check, P.available)), row.names = FALSE)

utils::write.csv(gamma_check,
                 file.path(PATHS$tables, "gamma_matrix_inventory.csv"),
                 row.names = FALSE)

