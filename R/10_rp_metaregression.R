# ======================================================================
# 10_rp_metaregression.R - z(Rp) ~ gamma_ij, AND SIGN CONCORDANCE
#
# Manuscript sections served:
#   Results > Bivariate model, final two paragraphs
#     "A complementary meta-regression of z(Rp) on gamma_ij ... slope =
#      0.113, 95% CI [0.020, 0.205], z = 2.39, P = 0.017; k = 535 trait
#      pairs from 38 studies"
#     "Counted by sign, the two were not independent: 55% of trait pairs
#      shared the sign ... (chi2_1 = 5.52, P = 0.019)"
#   Figure 2b (fitted line and quadrant percentages)
#
# Requires: df (01_load_clean.R)
# Produces: df_rp (shared with 11 and 12), rma_rp_z_full, rma_rp_z,
#           rma_rp_z_ss, mcmc_rp*, concordance
#
# EVERYTHING RUNS ON THE FISHER-z SCALE
#   Response is z(Rp) = atanh(Rp); sampling variance is 1/(n-3).
# ======================================================================


# -- 1. Analysis frame -------------------------------------------------

df_rp <- df %>%
  dplyr::filter(
    !is.na(Gamma.ij),
    !is.na(Gamma.ij.se), Gamma.ij.se > 0, is.finite(Gamma.ij.se),
    !is.na(Rp), abs(Rp) < 1,
    !is.na(Sample.size), Sample.size > 3,
    !is.na(Kingdom),
    !is.na(Type.of.fitness),
    !is.na(Type.of.selection),
    !is.na(Type.match),
    !is.na(Setting.field.lab.artificial),
    !is.na(ID),
    !is.na(Trait.2)
  )

df_rp <- df_rp %>%
  dplyr::mutate(
    dplyr::across(c(Kingdom, Type.of.fitness, Type.of.selection, ID,
                    Type.match, Setting.field.lab.artificial, group_id),
                  as.factor),
    Sample.size.c = log(Sample.size)   # log sample size; name kept for continuity
  )

df_rp$Rp_z   <- atanh(df_rp$Rp)
df_rp$Vz     <- 1 / (df_rp$Sample.size - 3)
df_rp$obs_id <- factor(seq_len(nrow(df_rp)))

# Conservative sampling variance for the sensitivity fit.
# --- OPEN QUESTION, see README "Known discrepancies" ---
# A design effect DE = 1 + (M - 1) * ICC inflates the VARIANCE by DE.
# With M = 2 and ICC = 0.5 that is Vz * 1.5. The original script used
# Vz * sqrt(1.5), which inflates the standard error by that factor
# instead. Set DESIGN_EFFECT_ON_VARIANCE to TRUE once the intended
# meaning is settled; FALSE reproduces the current manuscript.
DESIGN_EFFECT_ON_VARIANCE <- FALSE
de <- 1 + (2 - 1) * 0.5
df_rp$Vz_con <- if (DESIGN_EFFECT_ON_VARIANCE) df_rp$Vz * de else df_rp$Vz * sqrt(de)

cat("\n---- 10_rp_metaregression ----\n")
cat("Trait pairs:", nrow(df_rp), " (Results: 535)\n")
cat("Studies    :", nlevels(droplevels(df_rp$ID)), " (Results: 38)\n")
stopifnot(all(abs(df_rp$Rp) < 1), all(df_rp$Sample.size > 3))

RE3_rp <- list(~ 1 | ID, ~ 1 | group_id, ~ 1 | obs_id)


# ======================================================================
# 2. METAFOR (REML) - PRIMARY ENGINE
# ======================================================================

# 2a. Full model, all moderators. Results: "None of the other
# predictors, kingdom, type of selection, trait-type match, fitness
# measure, or experimental setting, explained additional variation".
rma_rp_z_full <- metafor::rma.mv(
  yi = Rp_z, V = Vz,
  mods = ~ Gamma.ij + Sample.size.c + Kingdom + Type.match +
           Type.of.fitness + Type.of.selection + Setting.field.lab.artificial,
  random = RE3_rp, data = df_rp, method = "REML")

cat("\n========== rma.mv: full ==========\n")
print(summary(rma_rp_z_full))

# 2b. Reduced model, gamma only. This is the slope quoted first in the
# Results (0.113, P = 0.017).
rma_rp_z <- metafor::rma.mv(
  yi = Rp_z, V = Vz, mods = ~ Gamma.ij,
  random = RE3_rp, data = df_rp, method = "REML")

cat("\n========== rma.mv: reduced, z(Rp) ~ gamma (REPORTED) ==========\n")
print(summary(rma_rp_z))

# 2c. With sample size retained. Results: "The effect of study sample
# size was also a significant negative term (P = 0.003) ... but the
# gamma_ij slope was robust to its inclusion (0.131, P = 0.005)."
rma_rp_z_ss <- metafor::rma.mv(
  yi = Rp_z, V = Vz, mods = ~ Gamma.ij + Sample.size.c,
  random = RE3_rp, data = df_rp, method = "REML")

cat("\n========== rma.mv: gamma + log(sample size) (REPORTED) ==========\n")
print(summary(rma_rp_z_ss))

# 2d. Sensitivity: conservative sampling variance.
# --- CHECK BEFORE SUBMISSION, see README ---
# Figure 2b draws its line and its slope annotation from THIS model, not
# from rma_rp_z. The two give a very similar slope, but the figure
# caption and the Results sentence should name the same model.
rma_rp_z_con <- metafor::rma.mv(
  yi = Rp_z, V = Vz_con, mods = ~ Gamma.ij + Sample.size.c,
  random = RE3_rp, data = df_rp, method = "REML")

cat("\n========== rma.mv: conservative variance (drives Figure 2b) ==========\n")
print(summary(rma_rp_z_con))

slope_compare <- data.frame(
  model = c("gamma only (Vz)", "gamma + log(n) (Vz)",
            "gamma + log(n) (Vz_con), used by Figure 2b"),
  slope = round(c(rma_rp_z$b["Gamma.ij", 1],
                  rma_rp_z_ss$b["Gamma.ij", 1],
                  rma_rp_z_con$b["Gamma.ij", 1]), 4),
  p     = signif(c(rma_rp_z$pval[2], rma_rp_z_ss$pval[2],
                   rma_rp_z_con$pval[2]), 3)
)
cat("\n---- gamma_ij slope across specifications ----\n")
print(slope_compare, row.names = FALSE)


# ======================================================================
# 3. MCMCglmm CROSS-CHECK
# ======================================================================
# Methods: "The same model specification was additionally fitted in a
# Bayesian framework using MCMCglmm as a cross-validation."
# Results: "with a Bayesian fit in close agreement (MCMCglmm slope =
# 0.115, pMCMC = 0.01)."
#
# Slow. Set RUN_MCMC to FALSE to skip when only the REML results are
# needed.

RUN_MCMC <- TRUE

if (RUN_MCMC) {
  prior_2G <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 0.002),
                            G2 = list(V = 1, nu = 0.002)))

  mcmc_args <- list(random = ~ ID + group_id, mev = df_rp$Vz, data = df_rp,
                    family = "gaussian", prior = prior_2G,
                    nitt = 1100000, burnin = 100000, thin = 1000,
                    verbose = FALSE)

  mcmc_rp <- do.call(MCMCglmm::MCMCglmm, c(list(fixed =
    Rp_z ~ Gamma.ij + Sample.size.c + Kingdom + Type.match +
           Type.of.fitness + Type.of.selection + Setting.field.lab.artificial),
    mcmc_args))
  cat("\n========== MCMCglmm: full ==========\n"); print(summary(mcmc_rp))

  mcmc_rp_red <- do.call(MCMCglmm::MCMCglmm,
                         c(list(fixed = Rp_z ~ Gamma.ij), mcmc_args))
  cat("\n========== MCMCglmm: reduced (REPORTED) ==========\n")
  print(summary(mcmc_rp_red))

  save_model(list(full = mcmc_rp, reduced = mcmc_rp_red), "10_mcmcglmm_rp")
} else {
  cat("\nRUN_MCMC is FALSE: skipping the MCMCglmm cross-check.\n")
}


# ======================================================================
# 4. SIGN CONCORDANCE
# ======================================================================
# Results: "55% of trait pairs shared the sign of the correlational
# gradient and the phenotypic correlation, against 45% that did not
# (chi2_1 = 5.52, P = 0.019)."
#
# This test runs on the same rows that Figure 2b plots: all trait pairs
# with a usable Rp and |gamma_ij| <= 1, excluding exact zeros in either
# quantity, which have no sign. It is a wider set than df_rp, which
# additionally requires every moderator to be present.

dB <- df %>%
  dplyr::filter(!is.na(Type.of.selection), !is.na(Rp), abs(Rp) < 1,
                !is.na(Gamma.ij), abs(Gamma.ij) <= 1)
dQ <- dB %>% dplyr::filter(Gamma.ij != 0, Rp != 0)

dQ$concordant <- sign(dQ$Gamma.ij) == sign(dQ$Rp)
n_conc <- sum(dQ$concordant); n_tot <- nrow(dQ)

cat("\n========== SIGN CONCORDANCE ==========\n")
cat(sprintf("Concordant: %d / %d (%.1f%%)\n", n_conc, n_tot,
            100 * n_conc / n_tot))

chi_test <- stats::chisq.test(c(n_conc, n_tot - n_conc), p = c(0.5, 0.5))
print(chi_test)

# One-sided exact test of the directional prediction
print(stats::binom.test(n_conc, n_tot, p = 0.5, alternative = "greater"))

# Study-clustered version. Trait pairs from one study are not
# independent, so this is the defensible test; the chi-square above is
# the one quoted, and the two should agree in direction.
if (requireNamespace("lme4", quietly = TRUE)) {
  m_conc <- lme4::glmer(concordant ~ 1 + (1 | ID), data = dQ, family = "binomial")
  print(summary(m_conc))
  cat("Clustering-adjusted P(concordant):",
      round(stats::plogis(lme4::fixef(m_conc)[["(Intercept)"]]), 3), "\n")
}

# Quadrant percentages, as annotated on Figure 2b
quadrants <- data.frame(
  quadrant = c("gamma > 0, Rp > 0", "gamma < 0, Rp > 0",
               "gamma > 0, Rp < 0", "gamma < 0, Rp < 0"),
  pct = round(100 * c(mean(dQ$Gamma.ij > 0 & dQ$Rp > 0),
                      mean(dQ$Gamma.ij < 0 & dQ$Rp > 0),
                      mean(dQ$Gamma.ij > 0 & dQ$Rp < 0),
                      mean(dQ$Gamma.ij < 0 & dQ$Rp < 0)), 1)
)
quadrants$type <- c("concordant", "discordant", "discordant", "concordant")
cat("\n---- Quadrant occupancy ----\n"); print(quadrants, row.names = FALSE)

concordance <- list(n_concordant = n_conc, n_total = n_tot,
                    pct = 100 * n_conc / n_tot,
                    chisq = chi_test, quadrants = quadrants)

utils::write.csv(quadrants,
                 file.path(PATHS$tables, "quadrant_occupancy.csv"),
                 row.names = FALSE)

save_model(list(full = rma_rp_z_full, reduced = rma_rp_z,
                with_n = rma_rp_z_ss, conservative = rma_rp_z_con,
                slope_compare = slope_compare,
                concordance = concordance),
           "10_rp_metaregression")
