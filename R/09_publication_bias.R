# ======================================================================
# 09_publication_bias.R - EGGER REGRESSIONS
#
# Manuscript sections served:
#   Methods > Heterogeneity and publication bias
#   Supplementary Figure S4 (funnel plots; drawn in
#     figures/figS4_funnel_plots.R from the models fitted here)
#
#
# Requires: df (01_load_clean.R)
# Produces: egger_gij_signed, egger_gij_abs, egger_gii_signed,
#           egger_gii_abs, egger_table

# ======================================================================


# -- 1. Prepare ---------------------------------------------------------

df_gij <- df[!is.na(df$Gamma.ij) & !is.na(df$Gamma.ij.se), ]
df_gii <- df[!is.na(df$Gamma.ii) & !is.na(df$Gamma.ii.se), ]

# Floor any zero standard errors at the smallest observed positive value;
# a zero SE would give an undefined sampling variance.
floor_se <- function(x) {
  if (any(x == 0, na.rm = TRUE)) x[x == 0] <- min(x[x > 0], na.rm = TRUE)
  x
}
df_gij$Gamma.ij.se <- floor_se(df_gij$Gamma.ij.se)
df_gii$Gamma.ii.se <- floor_se(df_gii$Gamma.ii.se)

df_gij$abs_gij <- abs(df_gij$Gamma.ij)
df_gii$abs_gii <- abs(df_gii$Gamma.ii)

df_gij$obs_id <- factor(seq_len(nrow(df_gij)))
df_gii$obs_id <- factor(seq_len(nrow(df_gii)))

for (v in c("ID", "group_id")) {
  df_gij[[v]] <- as.factor(df_gij[[v]])
  df_gii[[v]] <- as.factor(df_gii[[v]])
}

cat("\n---- 09_publication_bias ----\n")
cat("gamma_ij: k =", nrow(df_gij), "| studies:", nlevels(df_gij$ID),
    "| groups:", nlevels(df_gij$group_id), "\n")
cat("gamma_ii: k =", nrow(df_gii), "| studies:", nlevels(df_gii$ID),
    "| groups:", nlevels(df_gii$group_id), "\n")

RE3_egger <- list(~ 1 | ID, ~ 1 | group_id, ~ 1 | obs_id)


# -- 2. The four Egger regressions -------------------------------------

egger_gij_signed <- metafor::rma.mv(
  yi = Gamma.ij, V = Gamma.ij.se^2, mods = ~ Gamma.ij.se,
  random = RE3_egger, data = df_gij, method = "REML")

egger_gij_abs <- metafor::rma.mv(
  yi = abs_gij, V = Gamma.ij.se^2, mods = ~ Gamma.ij.se,
  random = RE3_egger, data = df_gij, method = "REML")

egger_gii_signed <- metafor::rma.mv(
  yi = Gamma.ii, V = Gamma.ii.se^2, mods = ~ Gamma.ii.se,
  random = RE3_egger, data = df_gii, method = "REML")

egger_gii_abs <- metafor::rma.mv(
  yi = abs_gii, V = Gamma.ii.se^2, mods = ~ Gamma.ii.se,
  random = RE3_egger, data = df_gii, method = "REML")

cat("\n========== Egger: signed gamma_ij (REPORTED) ==========\n")
print(summary(egger_gij_signed))
cat("\n========== Egger: signed gamma_ii (REPORTED) ==========\n")
print(summary(egger_gii_signed))
cat("\n========== Egger: |gamma_ij| (text only) ==========\n")
print(summary(egger_gij_abs))
cat("\n========== Egger: |gamma_ii| (text only) ==========\n")
print(summary(egger_gii_abs))


# -- 3. Summary table --------------------------------------------------

egger_row <- function(mod, label, k, n_studies) {
  data.frame(
    model     = label,
    k         = k,
    studies   = n_studies,
    slope     = round(mod$b[2, 1], 3),
    ci_lower  = round(mod$ci.lb[2], 3),
    ci_upper  = round(mod$ci.ub[2], 3),
    z         = round(mod$zval[2], 2),
    p         = signif(mod$pval[2], 3)
  )
}

egger_table <- rbind(
  egger_row(egger_gij_signed, "signed gamma_ij", nrow(df_gij), nlevels(df_gij$ID)),
  egger_row(egger_gii_signed, "signed gamma_ii", nrow(df_gii), nlevels(df_gii$ID)),
  egger_row(egger_gij_abs,    "|gamma_ij|",      nrow(df_gij), nlevels(df_gij$ID)),
  egger_row(egger_gii_abs,    "|gamma_ii|",      nrow(df_gii), nlevels(df_gii$ID))
)

cat("\n========== EGGER SUMMARY ==========\n")
print(egger_table, row.names = FALSE)

utils::write.csv(egger_table,
                 file.path(PATHS$tables, "egger_tests.csv"), row.names = FALSE)

save_model(list(gij_signed = egger_gij_signed, gij_abs = egger_gij_abs,
                gii_signed = egger_gii_signed, gii_abs = egger_gii_abs,
                data_gij = df_gij, data_gii = df_gii,
                table = egger_table),
           "09_publication_bias")
