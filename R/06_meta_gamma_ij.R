# ======================================================================
# 06_meta_gamma_ij.R - MAIN META-ANALYSIS OF gamma_ij
#
# Manuscript sections served:
#   Methods > Meta-analytic model
#   Methods > Heterogeneity (variance components)
#   Results > Strength and generality of correlational selection
#
# Requires: combined_df (05_build_combined.R)
#           rve_summary, variance_partition, fold_norm, save_model
#           (00_functions.R)
# Produces: model_signed, E_overall, level_table
#
# DESIGN
#   Every model is fitted on the SIGNED gradient yi with the diagonal
#   sampling variance vi (squared standard error from Masterdata).
#   No model uses |gamma_ij| as the response. The magnitude of
#   correlational selection is obtained from the signed model with the
#   folded-normal expectation (fold_norm), which uses the signed mean
#   and the summed random-effect variances.
#
#   Non-independence among effect sizes within a study is handled by the
#   study and group random effects and by cluster-robust (CR2) inference
#   with study as the cluster.
# ======================================================================

cat("\n---- 06_meta_gamma_ij ----\n")
cat("N effect sizes:", nrow(combined_df), "\n")
cat("N studies     :", dplyr::n_distinct(combined_df$study_id), "\n")
cat("N groups      :", dplyr::n_distinct(combined_df$group_id), "\n")

RE3 <- list(~ 1 | study_id, ~ 1 | group_id, ~ 1 | obs_id)


# ======================================================================
# 1. SIGNED INTERCEPT-ONLY MODEL
# ======================================================================
# The overall signed mean gamma_ij and the variance partition across
# the three random levels: study identity, within-study grouping, and
# the observation-level term that absorbs residual between-effect-size
# heterogeneity (Nakagawa and Santos 2012).

model_signed <- metafor::rma.mv(
  yi = yi, V = vi, random = RE3,
  data = combined_df, method = "REML")

cat("\n========== SIGNED INTERCEPT-ONLY ==========\n")
print(model_signed)
rve_summary(model_signed, combined_df$study_id, "signed intercept-only")

cat("\n-- Variance partition --\n")
variance_partition(model_signed,
                   labels = c("study", "group", "observation"))

cat("\n-- Heterogeneity --\n")
cat("Q =", round(model_signed$QE, 1),
    "on", model_signed$k - model_signed$p, "df,  P =",
    format.pval(model_signed$QEp, digits = 3), "\n")


# ======================================================================
# 2. OVERALL BIAS-CORRECTED MAGNITUDE
# ======================================================================
# E[|gamma_ij|] from the signed mean and the summed variance components
# of model_signed. This is a point estimate. Its bootstrap interval is
# computed in 08_bias_correction.R.

E_overall <- fold_norm(model_signed)

cat("\n-- Bias-corrected magnitude, all effect sizes --\n")
cat("signed mean     :", round(as.numeric(model_signed$beta[1]), 4), "\n")
cat("total tau       :", round(sqrt(sum(model_signed$sigma2)), 4), "\n")
cat("E[|gamma_ij|]   :", round(E_overall, 4), "\n")


# ======================================================================
# 3. BIAS-CORRECTED MAGNITUDE PER PREDICTOR LEVEL
# ======================================================================
# For each level of each categorical predictor, the signed model is
# fitted to that subset and the folded-normal expectation is taken from
# it. A predictor placed in a meta-regression on the signed response
# would test a shift in the mean SIGN of gamma_ij, not in its magnitude,
# so magnitude is compared between levels through separate fits.
#
# These are point estimates. Levels resting on few studies have poorly
# estimated variance components and are flagged. Type of selection is
# handled in 08_bias_correction.R, with bootstrap intervals.

level_magnitude <- function(var) {
  levs <- sort(unique(stats::na.omit(as.character(combined_df[[var]]))))
  do.call(rbind, lapply(levs, function(l) {
    d <- combined_df[!is.na(combined_df[[var]]) &
                       as.character(combined_df[[var]]) == l, ]
    fit <- tryCatch(
      metafor::rma.mv(yi, vi, random = RE3, data = d, method = "REML"),
      error = function(e) NULL)
    n_st <- dplyr::n_distinct(d$study_id)
    data.frame(
      predictor   = var,
      level       = l,
      k           = nrow(d),
      studies     = n_st,
      signed_mean = if (is.null(fit)) NA_real_ else round(as.numeric(fit$beta[1]), 4),
      tau_total   = if (is.null(fit)) NA_real_ else round(sqrt(sum(fit$sigma2)), 4),
      E_abs       = if (is.null(fit)) NA_real_ else round(fold_norm(fit), 4),
      note        = if (is.null(fit)) "model failed"
      else if (n_st < 10) "fewer than 10 studies" else "",
      stringsAsFactors = FALSE
    )
  }))
}

level_predictors <- c("Kingdom", "Setting", "Type.of.fitness", "Type.match")

level_table <- do.call(rbind, lapply(level_predictors, level_magnitude))

cat("\n========== BIAS-CORRECTED MAGNITUDE PER LEVEL ==========\n")
print(level_table, row.names = FALSE)


# ======================================================================
# 4. SAVE
# ======================================================================

utils::write.csv(level_table,
                 file.path(PATHS$tables, "magnitude_per_level.csv"),
                 row.names = FALSE)

save_model(list(model_signed = model_signed,
                E_overall    = E_overall,
                level_table  = level_table),
           "06_meta_gamma_ij")