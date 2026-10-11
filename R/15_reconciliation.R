# ======================================================================
# 15_reconciliation.R - EVERY REPORTED NUMBER AGAINST ITS COMPUTED VALUE
#
# This script produces nothing new. It walks the manuscript and the
# supplement, pulls the corresponding value out of the fitted objects,
# and prints a pass or fail against a tolerance. Run it last, every
# time, before touching the Overleaf document.
#
# Requires: everything from 01 to 14. Objects that are missing (because
#   a slow section was skipped) are reported as "not run" rather than
#   causing an error.
# Produces: output/tables/reconciliation.csv
#
# HOW TO READ THE OUTPUT
#   OK      computed value matches the manuscript within tolerance
#   CHECK   it does not; either the manuscript or the pipeline is stale
#   skipped the object was not in the workspace
#
# Five entries are expected to come back CHECK against the current
# Masterdata.xlsx. They are described in README "Known discrepancies";
# do not silently edit the manuscript to match, work out which side is
# wrong first.
# ======================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

.checks <- list()

check <- function(label, reported, computed, tol = 0.5, section = "") {
  val <- tryCatch(computed, error = function(e) NA_real_)
  status <- if (length(val) != 1 || is.na(val)) "skipped"
            else if (abs(val - reported) <= tol) "OK" else "CHECK"
  .checks[[length(.checks) + 1]] <<- data.frame(
    section  = section,
    quantity = label,
    reported = reported,
    computed = if (length(val) == 1 && !is.na(val)) round(val, 4) else NA,
    status   = status,
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

exists_obj <- function(x) exists(x, envir = .GlobalEnv, inherits = FALSE)

cat("\n")
cat("======================================================================\n")
cat(" RECONCILIATION: manuscript against pipeline\n")
cat("======================================================================\n")


# ---- Abstract and Results: sample --------------------------------------

check("total gradients entering analyses", 1649,
      sum(!is.na(df$Gamma.ii)) + sum(!is.na(df$Gamma.ij)), 0, "Results sample")
check("quadratic gradients (gamma_ii)", 639,
      sum(!is.na(df$Gamma.ii)), 0, "Results sample")
check("correlational gradients (gamma_ij)", 1010,
      sum(!is.na(df$Gamma.ij)), 0, "Results sample")
check("studies", 70, dplyr::n_distinct(df$ID), 0, "Results sample")
check("species", 63,
      dplyr::n_distinct(df$Latin.name[!is.na(df$Latin.name)]), 0,
      "Results sample")
check("taxonomic classes", 13,
      dplyr::n_distinct(df$Class[!is.na(df$Class) & df$Class != "N/A"]), 0,
      "Results sample")
check("orders", 38,
      dplyr::n_distinct(df$Order[!is.na(df$Order) & df$Order != "N/A"]), 0,
      "Results sample")

if (exists_obj("gamma_check")) {
  check("gamma-matrices", 213, nrow(gamma_check), 0, "Methods extraction")
  check("full-rank gamma-matrices", 106, sum(gamma_check$full_rank), 0,
        "Methods extraction")
}
if (exists_obj("tree")) {
  check("tips in the phylogeny (Figure 1)", 63, length(tree$tip.label), 0,
        "Figure 1")
}


# ---- Results: strength and generality ----------------------------------

if (exists_obj("E_nat")) {
  check("bias-corrected E[|gamma_ij|], natural", 0.087, E_nat, 0.005,
        "Results strength")
  check("bias-corrected E[|gamma_ij|], sexual", 0.234, E_sex, 0.02,
        "Results strength")
  check("ratio sexual / natural", 2.7, E_sex / E_nat, 0.2, "Results strength")
}
if (exists_obj("model_red3_noint")) {
  b <- as.numeric(model_red3_noint$beta)
  nm <- rownames(model_red3_noint$beta)
  check("raw |gamma_ij|, natural", 0.115,
        b[grep("Natural", nm)], 0.005, "Results strength")
  check("raw |gamma_ij|, sexual", 0.217,
        b[grep("Sexual", nm)], 0.01, "Results strength")
}
if (exists_obj("model_intercept_vi")) {
  check("QE, main gamma_ij model", 8917.8, model_intercept_vi$QE, 50,
        "Results strength")
  s2 <- model_intercept_vi$sigma2
  check("variance at study level, per cent", 42, 100 * s2[1] / sum(s2), 3,
        "Results strength")
  check("variance at group level, per cent", 37, 100 * s2[2] / sum(s2), 3,
        "Results strength")
  check("variance at observation level, per cent", 20,
        100 * s2[3] / sum(s2), 3, "Results strength")
}


# ---- Results: bivariate model -------------------------------------------

if (exists_obj("df_rp")) {
  check("trait pairs in the Rp analyses", 535, nrow(df_rp), 0,
        "Results bivariate")
  check("studies in the Rp analyses", 38,
        dplyr::n_distinct(df_rp$ID), 0, "Results bivariate")
}
if (exists_obj("fit_biv")) {
  check("signed mean gamma_ij (bivariate)", 0.013,
        as.numeric(fit_biv$beta)[1], 0.005, "Results bivariate")
  check("mean z(Rp) (bivariate)", 0.142,
        as.numeric(fit_biv$beta)[2], 0.01, "Results bivariate")
  check("REML rho_B", 0.78, fit_biv$rho, 0.05, "Results bivariate")
}
if (exists_obj("rho_summary")) {
  check("Bayesian rho_B (median)", 0.63, rho_summary$median, 0.05,
        "Results bivariate")
  check("posterior P(rho_B > 0)", 0.97, rho_summary$p_gt_zero, 0.03,
        "Results bivariate")
}
if (exists_obj("rma_rp_z")) {
  check("meta-regression slope of z(Rp) on gamma_ij", 0.113,
        as.numeric(rma_rp_z$beta[2]), 0.01, "Results bivariate")
  check("its P value", 0.017, rma_rp_z$pval[2], 0.01, "Results bivariate")
}
if (exists_obj("rma_rp_z_ss")) {
  check("slope with log sample size retained", 0.131,
        as.numeric(rma_rp_z_ss$beta[2]), 0.01, "Results bivariate")
}
if (exists_obj("concordance")) {
  check("concordant trait pairs, per cent", 55, concordance$pct, 1,
        "Results bivariate")
  check("chi-square for sign concordance", 5.52,
        as.numeric(concordance$chisq$statistic), 0.2, "Results bivariate")
  check("its P value", 0.019, as.numeric(concordance$chisq$p.value), 0.005,
        "Results bivariate")
}


# ---- Results: angles ----------------------------------------------------

if (exists_obj("all_angles")) {
  check("gamma-P comparisons", 75, nrow(all_angles), 0, "Results angles")
  check("studies contributing comparisons", 26,
        dplyr::n_distinct(all_angles$StudyID), 0, "Results angles")
  check("comparisons with k >= 2", 24,
        sum(all_angles$Parameters >= 2), 0, "Results angles")
  d <- with(all_angles[all_angles$Parameters >= 2, ], Angle_max - Angle)
  check("mean gain at k >= 2, degrees (Figure 3a)", 45.7,
        mean(d, na.rm = TRUE), 2, "Figure 3a")
}
if (exists_obj("angle_full")) {
  check("angle intercept, degrees", 48.5, as.numeric(angle_full$beta[1]), 1,
        "Results angles")
  check("angle slope, degrees per unit eigenvalue", -20.2,
        as.numeric(angle_full$beta[2]), 1.5, "Results angles")
  check("sigma2 study (angle model)", 203.3, angle_full$sigma2[1], 15,
        "Results angles")
  check("sigma2 observation (angle model)", 281.2, angle_full$sigma2[2], 15,
        "Results angles")
  check("QE (angle model)", 5185.7, angle_full$QE, 100, "Results angles")
}
if (exists_obj("angle_cr2")) {
  check("CR2 Satterthwaite df for the angle slope", 3.02,
        angle_cr2$df_Satt[2], 0.5, "Results angles")
  check("CR2 t for the angle slope", -4.67, angle_cr2$tstat[2], 0.5,
        "Results angles")
}
if (exists_obj("sim_full")) {
  check("similarity intercept, logit scale", -0.81,
        as.numeric(sim_full$beta[1]), 0.1, "Results angles")
  check("similarity slope, logit scale", 0.83,
        as.numeric(sim_full$beta[2]), 0.15, "Results angles")
  check("QE (similarity model)", 1745.4, sim_full$QE, 100, "Results angles")
}


# ---- Supplement: Egger tests --------------------------------------------

if (exists_obj("egger_gij_signed")) {
  check("Egger slope, signed gamma_ij", 0.188, egger_gij_signed$b[2, 1], 0.02,
        "Figure S4")
  check("its P value", 0.078, egger_gij_signed$pval[2], 0.02, "Figure S4")
  check("Egger slope, signed gamma_ii", 0.132, egger_gii_signed$b[2, 1], 0.02,
        "Figure S4")
  check("its P value", 0.319, egger_gii_signed$pval[2], 0.05, "Figure S4")
}


# ---- Output --------------------------------------------------------------

recon <- do.call(rbind, .checks)

cat("\n")
print(recon, row.names = FALSE)

n_ok      <- sum(recon$status == "OK")
n_check   <- sum(recon$status == "CHECK")
n_skipped <- sum(recon$status == "skipped")

cat("\n----------------------------------------------------------------------\n")
cat(sprintf(" %d OK   |   %d CHECK   |   %d skipped\n",
            n_ok, n_check, n_skipped))
cat("----------------------------------------------------------------------\n")

if (n_check > 0) {
  cat("\nEntries needing attention:\n")
  print(recon[recon$status == "CHECK", ], row.names = FALSE)
  cat("\nCompare against README section \"Known discrepancies\" before\n")
  cat("changing anything in the manuscript.\n")
}

utils::write.csv(recon, file.path(PATHS$tables, "reconciliation.csv"),
                 row.names = FALSE)
cat("\nWritten to output/tables/reconciliation.csv\n")
