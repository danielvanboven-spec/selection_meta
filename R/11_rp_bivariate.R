# ======================================================================
# 11_rp_bivariate.R - BIVARIATE RANDOM-EFFECTS META-ANALYSIS (metafor)
#
# Manuscript sections served:
#   Methods > Bivariate model of phenotypic correlations and
#             correlational selection gradients
#   Results > Bivariate model
#     signed mean gamma_ij = 0.013, 95% CI [-0.011, 0.038], P = 0.29
#     mean z(Rp)           = 0.142, 95% CI [ 0.022, 0.262], P = 0.021
#     REML rho_B           = 0.78, weakly identified
#
# Requires: df_rp (10_rp_metaregression.R)
# Produces: fit_biv, rho_sensitivity
#
# WHY A BIVARIATE MODEL RATHER THAN A META-REGRESSION
#   A meta-regression of Rp on gamma_ij treats gamma_ij as a fixed,
#   error-free predictor, so its sampling error is discarded and the
#   relationship is attenuated by regression dilution. It also assumes
#   an independent predictor, which fails because gamma_ij and Rp are
#   estimated from the same individuals. The bivariate model treats both
#   as outcomes, each carrying its own sampling variance, and estimates
#   the between-study correlation of their true effects, rho_B, directly.
#
# MODEL
#   y_kj = mu_j + u_kj + e_kj + eps_kj
#   with (u_k,gamma, u_k,z) ~ N(0, [tau_g^2, rho_B tau_g tau_z;
#                                   rho_B tau_g tau_z, tau_z^2])
#   ~ 0 + outcome              gives an outcome-specific grand mean
#   ~ outcome | study_id, UN   estimates tau_g^2, tau_z^2 and rho_B
#
# V
#   One 2x2 block per trait pair. Diagonal entries are the sampling
#   variances of gamma_hat_ij and z(Rp_hat); the off-diagonal is
#   r_w sqrt(Var(gamma) Var(z)). r_w cannot be derived analytically
#   without the raw trait covariances, so it is set to R_WITHIN (0.5)
#   for the primary fit and varied over R_WITHIN_SENS in section 5.
#
# THE REML rho_B SITS NEAR THE PARAMETER BOUNDARY
#   Its profile interval runs to about [0.00, 0.98], so the Results lead
#   with the Bayesian estimate from 12_rp_bivariate_brms.R and report
#   this one as agreeing in direction. Keep that ordering.
# ======================================================================


# -- 1. Data -----------------------------------------------------------

dat          <- df_rp
dat$gamma_vi <- dat$Gamma.ij.se^2
dat$zRp      <- dat$Rp_z
dat$zRp_var  <- dat$Vz

cat("\n---- 11_rp_bivariate ----\n")
cat("Trait pairs:", nrow(dat), "\n")
cat("Studies    :", dplyr::n_distinct(dat$ID), "\n")


# -- 2. Long format: two rows per trait pair ---------------------------

long <- dat %>%
  dplyr::mutate(row_k = dplyr::row_number()) %>%
  dplyr::transmute(
    row_k,
    study_id = as.character(ID),
    group_id = as.character(group_id),
    pair_id  = row_k,
    gamma_yi = Gamma.ij, gamma_vi = gamma_vi,
    zRp_yi   = zRp,      zRp_vi   = zRp_var) %>%
  tidyr::pivot_longer(c(gamma_yi, gamma_vi, zRp_yi, zRp_vi),
                      names_to = c("outcome", ".value"),
                      names_sep = "_") %>%
  dplyr::arrange(row_k, outcome) %>%
  dplyr::mutate(outcome = factor(outcome, levels = c("gamma", "zRp")))


# -- 3. Bivariate sampling-variance matrix -----------------------------

build_V_biv <- function(v_g, v_z, r_within) {
  k <- length(v_g)
  g_idx   <- seq(1, 2 * k, by = 2)
  zRp_idx <- seq(2, 2 * k, by = 2)
  V <- matrix(0, 2 * k, 2 * k)
  V[cbind(g_idx,   g_idx)]   <- v_g
  V[cbind(zRp_idx, zRp_idx)] <- v_z
  cov_gz <- r_within * sqrt(v_g * v_z)
  V[cbind(g_idx,   zRp_idx)] <- cov_gz
  V[cbind(zRp_idx, g_idx)]   <- cov_gz
  V
}

v_g <- dat$gamma_vi
v_z <- dat$zRp_var
V_biv <- build_V_biv(v_g, v_z, R_WITHIN)

stopifnot(isSymmetric(V_biv), !anyNA(V_biv),
          nrow(V_biv) == nrow(long))


# -- 4. Fit ------------------------------------------------------------

fit_biv <- metafor::rma.mv(
  yi     = yi,
  V      = V_biv,
  mods   = ~ 0 + outcome,
  random = list(~ outcome | study_id, ~ 1 | pair_id),
  struct = "UN",
  data   = long,
  method = "REML",
  sparse = TRUE,
  control = list(optimizer = "nlminb", rel.tol = 1e-8)
)

cat("\n========== BIVARIATE MODEL (r_within =", R_WITHIN, ") ==========\n")
print(summary(fit_biv))

cat("\nrho_B      =", round(fit_biv$rho, 3), "\n")
cat("tau_gamma  =", round(sqrt(fit_biv$tau2[1]), 3),
    "  tau_zRp =", round(sqrt(fit_biv$tau2[2]), 3), "\n")

# Outcome-specific grand means, as quoted in the Results
biv_means <- data.frame(
  outcome  = c("gamma_ij", "z(Rp)"),
  estimate = round(as.numeric(fit_biv$beta), 3),
  ci_lower = round(fit_biv$ci.lb, 3),
  ci_upper = round(fit_biv$ci.ub, 3),
  p        = signif(fit_biv$pval, 3)
)
cat("\n---- Outcome means ----\n"); print(biv_means, row.names = FALSE)

# Profile-likelihood intervals. These often run to the boundary for
# rho_B; that is the weak identification noted in the Results, not a
# convergence failure.
ci_biv <- tryCatch(stats::confint(fit_biv),
                   error = function(e) {
                     cat("confint failed:", conditionMessage(e), "\n"); NULL })
print(ci_biv)


# -- 5. Sensitivity to the within-pair sampling correlation ------------
# Methods: "I set r_w = 0.5 for the primary analysis and conducted
# sensitivity analyses over r_w in {0, 0.3, 0.5, 0.7, 0.9}".

cat("\n---- Sensitivity: r_within ----\n")
rho_sensitivity <- do.call(rbind, lapply(R_WITHIN_SENS, function(rw) {
  f <- tryCatch(
    metafor::rma.mv(yi, build_V_biv(v_g, v_z, rw), mods = ~ 0 + outcome,
                    random = list(~ outcome | study_id, ~ 1 | pair_id),
                    struct = "UN", data = long, method = "REML", sparse = TRUE,
                    control = list(optimizer = "nlminb", rel.tol = 1e-8)),
    error = function(e) NULL)
  data.frame(r_within = rw,
             rho_B = if (is.null(f)) NA_real_ else round(f$rho, 3))
}))
print(rho_sensitivity, row.names = FALSE)

utils::write.csv(rho_sensitivity,
                 file.path(PATHS$tables, "rho_B_sensitivity.csv"),
                 row.names = FALSE)

save_model(list(fit = fit_biv, means = biv_means,
                sensitivity = rho_sensitivity, ci = ci_biv),
           "11_rp_bivariate")
