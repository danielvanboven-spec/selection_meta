# ======================================================================
# 08_bias_correction.R - FOLDED-NORMAL BIAS CORRECTION
#
# Manuscript sections served:
#   Methods > Meta-analytic model, final two paragraphs
#   Results > Strength and generality of correlational selection
#   Abstract (E[|gamma_ij|] ~ 0.087)
#
# The mean of the absolute gradients overstates the magnitude of
# correlational selection, because |gamma_hat| carries sampling error
# whenever the underlying signed gradient is near zero. The magnitude is
# therefore taken from the SIGNED model and corrected with the
# folded-normal expectation (Morrissey 2016):
#
#   E[|gamma|] = tau sqrt(2/pi) exp(-mu^2 / 2 tau^2)
#                + mu (2 Phi(mu / tau) - 1)
#
# with mu the signed meta-analytic mean and tau^2 the sum of the
# random-effect variance components of the signed model. The correction
# is computed separately for natural and sexual selection, the one
# predictor retained by the backward elimination in 06.
#
# Numbers this script must produce:
#   natural E[|gamma_ij|] = 0.087, bootstrap 95% CI [0.036, 0.143]
#   sexual  E[|gamma_ij|] = 0.234, bootstrap 95% CI [0.063, 0.332]
#   ratio                 = 2.69
#
# Requires: combined_df (05_build_combined.R)
# Produces: fit_nat, fit_sex, E_nat, E_sex, boot_nat, boot_sex
#
# CAUTION ON THE SEXUAL ESTIMATE
#   The sexual subset rests on about ten study clusters. The study-level
#   variance collapses to the boundary there, so the two-level random
#   structure is used and the estimate is reported as a caution rather
#   than an established value. Keep that framing in the manuscript.
#
# RUNTIME
#   Two cluster bootstraps at B = N_BOOT. Expect several minutes.
# ======================================================================

RE3 <- list(~ 1 | study_id, ~ 1 | group_id, ~ 1 | obs_id)
RE2 <- list(~ 1 | study_id, ~ 1 | obs_id)


# -- 1. Signed model, overall ------------------------------------------
# Parallel to model_red3 in 06 but on the signed response.

model_red3_signed <- metafor::rma.mv(
  yi = yi, V = vi, random = RE3,
  data = combined_df, method = "REML")

cat("\n---- 08_bias_correction ----\n")
cat("\n========== SIGNED gamma_ij, overall ==========\n")
print(summary(model_red3_signed))

# The signed mean sits near zero, which is the whole reason the absolute
# mean needs correcting: positive and negative gradients are about
# equally common.
cat("\nSigned mean gamma_ij:", round(as.numeric(model_red3_signed$beta[1]), 4),
    " P =", format.pval(model_red3_signed$pval[1], digits = 3), "\n")


# -- 2. Split by type of selection -------------------------------------

nat <- combined_df[combined_df$Type.of.selection == "Natural selection", ]
sex <- combined_df[combined_df$Type.of.selection == "Sexual selection", ]

cat("\nnatural: k =", nrow(nat), "from", dplyr::n_distinct(nat$study_id), "studies\n")
cat("sexual : k =", nrow(sex), "from", dplyr::n_distinct(sex$study_id), "studies\n")

fit_nat <- metafor::rma.mv(yi, vi, random = RE3, data = nat, method = "REML")
fit_sex <- metafor::rma.mv(yi, vi, random = RE2, data = sex, method = "REML")

cat("\n========== SIGNED gamma_ij: natural ==========\n"); print(summary(fit_nat))
cat("\n========== SIGNED gamma_ij: sexual ==========\n");  print(summary(fit_sex))


# -- 3. Point estimates ------------------------------------------------

E_nat <- fold_norm(fit_nat)
E_sex <- fold_norm(fit_sex)

cat("\n---- Bias-corrected magnitudes ----\n")
cat("natural E[|gamma_ij|] :", round(E_nat, 3), "\n")
cat("sexual  E[|gamma_ij|] :", round(E_sex, 3), "\n")
cat("ratio (sexual/natural):", round(E_sex / E_nat, 2), "\n")


# -- 4. Cluster bootstrap confidence intervals -------------------------
# Methods: "Confidence intervals were obtained by a non-parametric
# cluster bootstrap, resampling studies with replacement, refitting the
# signed model and recomputing E[|gamma_ij|], with 95 per cent
# percentile intervals across B = 1,000 replicates."
#
# The original script ran B = 100 while the Methods claimed 1,000. It is
# now driven by N_BOOT in 00_setup.R so the two cannot diverge again.

cat("\nBootstrapping with B =", N_BOOT, "...\n")

boot_nat <- boot_fold(nat, B = N_BOOT, formula_random = RE3)
boot_sex <- boot_fold(sex, B = N_BOOT, formula_random = RE2)

ci_nat <- stats::quantile(boot_nat, c(0.025, 0.975), na.rm = TRUE)
ci_sex <- stats::quantile(boot_sex, c(0.025, 0.975), na.rm = TRUE)

n_fail_nat <- sum(is.na(boot_nat)); n_fail_sex <- sum(is.na(boot_sex))
if (n_fail_nat + n_fail_sex > 0)
  cat("Bootstrap refits that failed to converge: natural", n_fail_nat,
      "| sexual", n_fail_sex, "\n")

boot_ratio <- boot_sex / boot_nat
ci_ratio   <- stats::quantile(boot_ratio, c(0.025, 0.975), na.rm = TRUE)


# -- 5. Report ---------------------------------------------------------
# Values are printed from the objects, never hard-coded, so a rerun that
# shifts an estimate cannot leave a stale number on screen.

bias_table <- data.frame(
  selection = c("Natural", "Sexual"),
  E_abs     = round(c(E_nat, E_sex), 3),
  ci_lower  = round(c(ci_nat[1], ci_sex[1]), 3),
  ci_upper  = round(c(ci_nat[2], ci_sex[2]), 3),
  k         = c(nrow(nat), nrow(sex)),
  studies   = c(dplyr::n_distinct(nat$study_id),
                dplyr::n_distinct(sex$study_id))
)

cat("\n========== BIAS-CORRECTED E[|gamma_ij|] ==========\n")
print(bias_table, row.names = FALSE)
cat("\nratio =", round(E_sex / E_nat, 2),
    " 95% CI [", round(ci_ratio[1], 2), ",", round(ci_ratio[2], 2), "]\n")

utils::write.csv(bias_table,
                 file.path(PATHS$tables, "bias_corrected_magnitudes.csv"),
                 row.names = FALSE)

save_model(list(signed_overall = model_red3_signed,
                natural = fit_nat, sexual = fit_sex,
                E_nat = E_nat, E_sex = E_sex,
                boot_nat = boot_nat, boot_sex = boot_sex,
                table = bias_table),
           "08_bias_correction")
