# ======================================================================
# 07_gamma_ii.R - QUADRATIC GRADIENTS (gamma_ii)
#
# Manuscript sections served:
#   Methods > Study search and selection
#     "Quadratic gradients (gamma_ii), which are estimated alongside
#      correlational gradients within the same model, were extracted
#      where reported."
#   Results > Sample and study characteristics (the 639 count)
#   Discussion (quadratic selection as an internal benchmark for the
#     magnitude of correlational selection)
#
# Requires: df (01_load_clean.R)
# Produces: combined_df_ii, model_quad_signed, fit_nat_q, fit_sex_q
#
# gamma_ii rows are the rows with no second trait. Sampling variance is
# the squared reported standard error; V_combined does not apply here,
# because the block-diagonal covariances were built for the gamma_ij
# analysis only.
# ======================================================================


# -- 1. Build the quadratic analysis frame -----------------------------

quad_pub <- df %>%
  dplyr::filter(is.na(Trait.2), !is.na(Gamma.ii),
                !is.na(Gamma.ii.se), Gamma.ii.se > 0)

combined_df_ii <- data.frame(
  yi         = quad_pub$Gamma.ii,
  vi         = quad_pub$Gamma.ii.se^2,
  study_id   = as.character(quad_pub$ID),
  Latin.name = quad_pub$Latin.name,
  trait1     = quad_pub$Trait.1,
  type       = "quadratic",
  stringsAsFactors = FALSE
)

combined_df_ii$group_id <- paste0(combined_df_ii$study_id, "_",
                                  combined_df_ii$Latin.name, "_",
                                  quad_pub$Grouping, "_",
                                  quad_pub$Fitness.measurement)
combined_df_ii$obs_id <- as.character(seq_len(nrow(combined_df_ii)))
combined_df_ii$abs_yi <- abs(combined_df_ii$yi)


# -- 2. Attach moderators ----------------------------------------------

mod_lookup_ii <- df %>%
  dplyr::filter(!is.na(Gamma.ii)) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(
    Kingdom           = dplyr::first(Kingdom),
    Class             = dplyr::first(Class),
    Order             = dplyr::first(Order),
    Type.of.selection = dplyr::first(Type.of.selection),
    Setting           = dplyr::first(Setting.field.lab.artificial),
    Type.of.fitness   = dplyr::first(Type.of.fitness),
    Sample.size       = dplyr::first(Sample.size),
    No..of.parameters = dplyr::first(No..of.parameters),
    .groups = "drop"
  )

combined_df_ii$study_id_num <- as.numeric(gsub("_.*", "", combined_df_ii$study_id))
combined_df_ii <- merge(combined_df_ii, mod_lookup_ii,
                        by.x = "study_id_num", by.y = "ID",
                        all.x = TRUE, sort = FALSE)
combined_df_ii <- combined_df_ii[order(as.numeric(combined_df_ii$obs_id)), ]
rownames(combined_df_ii) <- NULL

cat("\n---- 07_gamma_ii ----\n")
cat("Quadratic gamma_ii rows:", nrow(combined_df_ii), " (Results: 639)\n")
cat("Studies                :", dplyr::n_distinct(combined_df_ii$study_id), "\n")


# -- 3. Signed model, overall ------------------------------------------

RE3 <- list(~ 1 | study_id, ~ 1 | group_id, ~ 1 | obs_id)
RE2 <- list(~ 1 | study_id, ~ 1 | obs_id)

model_quad_signed <- metafor::rma.mv(
  yi, vi, random = RE3, data = combined_df_ii, method = "REML")

cat("\n========== SIGNED gamma_ii ==========\n")
print(summary(model_quad_signed))
variance_partition(model_quad_signed, c("study", "group", "observation"))


# -- 4. Split by type of selection -------------------------------------
# The sexual-selection subset uses the two-level random structure. The
# group-level variance collapses to the boundary there for the same
# reason as in the gamma_ij models: within a group, signed gradients
# average to approximately zero, so there is no between-group variance
# left for the term to capture. This is a boundary collapse, not a
# failure of the model.

nat_q <- combined_df_ii[combined_df_ii$Type.of.selection == "Natural selection", ]
sex_q <- combined_df_ii[combined_df_ii$Type.of.selection == "Sexual selection", ]

fit_nat_q <- metafor::rma.mv(yi, vi, random = RE3, data = nat_q, method = "REML")
fit_sex_q <- metafor::rma.mv(yi, vi, random = RE2, data = sex_q, method = "REML")

cat("\n========== SIGNED gamma_ii: natural ==========\n"); print(summary(fit_nat_q))
cat("\n========== SIGNED gamma_ii: sexual ==========\n");  print(summary(fit_sex_q))


# -- 5. Bias-corrected magnitude ---------------------------------------
# Same folded-normal correction as for gamma_ij (08_bias_correction.R),
# giving the internal quadratic benchmark used in the Discussion.

cat("\n---- Folded-normal E[|gamma_ii|] ----\n")
cat("natural:", round(fold_norm(fit_nat_q), 3), "\n")
cat("sexual :", round(fold_norm(fit_sex_q), 3), "\n")

save_model(list(signed = model_quad_signed,
                natural = fit_nat_q, sexual = fit_sex_q,
                data = combined_df_ii),
           "07_gamma_ii")
