# ======================================================================
# 12_rp_bivariate_brms.R - BAYESIAN BIVARIATE META-ANALYSIS
#
# Manuscript sections served:
#   Results > Bivariate model
#     "A Bayesian fit, which estimates the within-pair residual
#      correlation directly, revealed that rho_B = 0.63 (95% CrI
#      [-0.07, 0.94]); almost the entire posterior was positive
#      (posterior probability rho_B > 0 = 0.97)."
#   Abstract (rho_B ~ 0.63)
#
# Requires: df_rp (10_rp_metaregression.R)
# Produces: fit_brms, rho_post, rho_summary
#
# WHY THIS FIT LEADS AND THE REML FIT FOLLOWS
#   11_rp_bivariate.R must be told the within-pair sampling correlation
#   r_w, because metafor cannot estimate it. Here set_rescor(TRUE)
#   estimates that residual correlation from the data, so rho_B is not
#   conditional on an assumed r_w. That is why the Results quote this
#   estimate first and report the REML one as agreeing in direction.
#
# STRUCTURE
#   (1 | p | study_id) makes the study-level random effects correlated
#   across the two outcomes, which is the brms equivalent of
#   struct = "UN" in metafor. That correlation is rho_B.
#
# DATA
#   Built from df_rp, the same rows as 10 and 11. An earlier version
#   drew on combined_df instead, which lost roughly ten studies' worth
#   of Rp in the merge, most of them sexual selection.
#
# RUNTIME
#   Four chains, 4,000 iterations each. Expect several minutes on four
#   cores. Set RUN_BRMS to FALSE to skip.
# ======================================================================

RUN_BRMS <- TRUE

if (!RUN_BRMS) {
  cat("\nRUN_BRMS is FALSE: skipping the Bayesian bivariate fit.\n")
} else {


# -- 1. Wide-format bivariate data -------------------------------------
# Each row carries both outcomes plus a known sampling SE for each.

dat_brms <- df_rp %>%
  dplyr::transmute(
    study_id          = factor(ID),
    group_id          = factor(group_id),
    Kingdom           = droplevels(factor(Kingdom)),
    Type.of.selection = droplevels(factor(Type.of.selection)),
    Type.match        = droplevels(factor(Type.match)),
    yi_gamma          = Gamma.ij,
    gamma_se          = Gamma.ij.se,   # reported sampling SE of gamma_ij
    yi_zRp            = Rp_z,          # atanh(Rp)
    zRp_se            = sqrt(Vz)       # sqrt(1/(n-3)), the Fisher-z SE
  )

cat("\n---- 12_rp_bivariate_brms ----\n")
cat("Trait pairs:", nrow(dat_brms), "\n")
cat("Studies    :", nlevels(droplevels(dat_brms$study_id)), "\n")


# -- 2. Model ----------------------------------------------------------

bf_gamma <- brms::bf(
  yi_gamma | se(gamma_se, sigma = TRUE) ~ 1 +
    (1 | p | study_id) + (1 | group_id))

bf_zRp <- brms::bf(
  yi_zRp | se(zRp_se, sigma = TRUE) ~ 1 +
    (1 | p | study_id) + (1 | group_id))

priors <- c(
  brms::prior(normal(0, 1),       class = "Intercept", resp = "yigamma"),
  brms::prior(normal(0, 1),       class = "Intercept", resp = "yizRp"),
  brms::prior(student_t(3, 0, 1), class = "sd",        resp = "yigamma"),
  brms::prior(student_t(3, 0, 1), class = "sd",        resp = "yizRp"),
  brms::prior(lkj(2),             class = "cor")
)

fit_brms <- brms::brm(
  formula = bf_gamma + bf_zRp + brms::set_rescor(TRUE),
  data    = dat_brms,
  prior   = priors,
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 10),
  seed    = SEED
)

cat("\n========== BAYESIAN BIVARIATE MODEL ==========\n")
print(summary(fit_brms))


# -- 3. Convergence diagnostics ----------------------------------------
# Check these before reporting anything. Rhat should be at or below
# 1.01 and the bulk effective sample size comfortably above 400.

post_summary <- posterior::summarise_draws(
  posterior::as_draws_df(fit_brms),
  "mean", "sd", "rhat", "ess_bulk", "ess_tail")
bad <- post_summary[!is.na(post_summary$rhat) & post_summary$rhat > 1.01, ]
if (nrow(bad) > 0) {
  cat("\nWARNING: parameters with Rhat > 1.01:\n"); print(bad)
} else {
  cat("\nAll Rhat <= 1.01.\n")
}


# -- 4. rho_B from the posterior ---------------------------------------

post <- posterior::as_draws_df(fit_brms)

cor_pars <- grep("^cor_study_id", names(post), value = TRUE)
cat("\nRandom-effect correlation parameters found:\n"); print(cor_pars)
stopifnot(length(cor_pars) >= 1)

rho_post <- post[[cor_pars[1]]]

rho_summary <- data.frame(
  parameter = "rho_B",
  median    = round(stats::median(rho_post), 3),
  crI_lower = round(stats::quantile(rho_post, 0.025), 3),
  crI_upper = round(stats::quantile(rho_post, 0.975), 3),
  p_gt_zero = round(mean(rho_post > 0), 3),
  row.names = NULL
)

cat("\n========== rho_B (REPORTED) ==========\n")
print(rho_summary, row.names = FALSE)


# -- 5. Residual correlation -------------------------------------------
# This is the within-pair sampling correlation that 11_rp_bivariate.R
# has to assume. Compare it with R_WITHIN in 00_setup.R: if the two
# disagree badly, the sensitivity table in 11 is the place to look.

rescor_pars <- grep("^rescor", names(post), value = TRUE)
if (length(rescor_pars) > 0) {
  rescor_post <- post[[rescor_pars[1]]]
  cat("\n---- Residual (within-pair) correlation ----\n")
  cat("median :", round(stats::median(rescor_post), 3), "\n")
  cat("95% CrI:", paste(round(stats::quantile(rescor_post, c(0.025, 0.975)), 3),
                        collapse = " to "), "\n")
  cat("r_within assumed in 11_rp_bivariate.R:", R_WITHIN, "\n")
}

utils::write.csv(rho_summary,
                 file.path(PATHS$tables, "rho_B_bayesian.csv"),
                 row.names = FALSE)

save_model(fit_brms, "12_fit_bivariate_brms")
save_model(list(rho_post = rho_post, summary = rho_summary),
           "12_rho_B_posterior")

}  # end RUN_BRMS
