# ======================================================================
# 05_build_combined.R - ASSEMBLE combined_df
#
# The analysis frame for every gamma_ij model: one row per correlational
# gradient in Masterdata, with its sampling variance (squared standard
# error) and the study-level moderators.
#
# Manuscript sections served:
#   Methods > Meta-analytic model
#
# Requires: df (01_load_clean.R)
# Produces: combined_df
#
# All values come from Masterdata. The raw-data refits of
# 04_lande_arnold_raw.R and the block-diagonal sampling matrix built
# from them are NOT used here. They enter only the sensitivity analysis
# in 16_sensitivity_raw_refits.R, so 04 does not have to be run for any
# reported result.
# ======================================================================


# -- 1. Correlational gradients from Masterdata ------------------------

gij <- df %>%
  dplyr::filter(!is.na(Gamma.ij), !is.na(Gamma.ij.se), Gamma.ij.se > 0)

combined_df <- data.frame(
  yi         = gij$Gamma.ij,
  vi         = gij$Gamma.ij.se^2,
  study_id   = as.character(gij$ID),
  Latin.name = gij$Latin.name,
  trait1     = gij$Trait.1,
  trait2     = gij$Trait.2,
  type       = "interaction",
  Rp         = gij$Rp,
  stringsAsFactors = FALSE
)

# group_id identifies the stratum (species x grouping x fitness measure)
# an estimate came from within its study
combined_df$group_id <- paste0(combined_df$study_id, "_",
                               combined_df$Latin.name, "_",
                               gij$Grouping, "_",
                               gij$Fitness.measurement)

combined_df$obs_id <- as.character(seq_len(nrow(combined_df)))
combined_df$abs_yi <- abs(combined_df$yi)


# -- 2. Attach moderators ----------------------------------------------
# Moderators are taken per study (first row of the study in Masterdata).

combined_df$study_id_num <- as.numeric(combined_df$study_id)

moderator_lookup <- df %>%
  dplyr::filter(!is.na(Gamma.ij)) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    Kingdom           = dplyr::first(Kingdom),
    Class             = dplyr::first(Class),
    Order             = dplyr::first(Order),
    Type.of.selection = dplyr::first(Type.of.selection),
    Setting           = dplyr::first(Setting.field.lab.artificial),
    Type.of.fitness   = dplyr::first(Type.of.fitness),
    Type.match        = dplyr::first(Type.match),
    Sample.size       = dplyr::first(Sample.size),
    No..of.parameters = dplyr::first(No..of.parameters),
    .groups = "drop"
  )

combined_df <- merge(combined_df, moderator_lookup,
                     by.x = "study_id_num", by.y = "ID",
                     all.x = TRUE, sort = FALSE)

# merge() can reorder rows; restore the Masterdata order
combined_df <- combined_df[order(as.numeric(combined_df$obs_id)), ]
rownames(combined_df) <- NULL

stopifnot(nrow(combined_df) == nrow(gij),
          all(as.integer(combined_df$obs_id) == seq_len(nrow(combined_df))))


# -- 3. Report ---------------------------------------------------------

cat("\n---- 05_build_combined ----\n")
cat("gamma_ij effect sizes   :", nrow(combined_df), "\n")
cat("Unique studies          :", dplyr::n_distinct(combined_df$study_id), "\n")
cat("Unique groups           :", dplyr::n_distinct(combined_df$group_id), "\n")
cat("Missing Latin.name      :", sum(is.na(combined_df$Latin.name)), "\n")

cat("\nNAs per moderator:\n")
print(sapply(combined_df[c("Kingdom", "Setting", "Type.of.selection",
                           "Type.match", "Type.of.fitness", "Sample.size")],
             function(x) sum(is.na(x))))