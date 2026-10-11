# ======================================================================
# 14_alignment_models.R - META-ANALYSIS OF SUBSPACE ALIGNMENT
#
# Manuscript sections served:
#   Methods > Subspace comparison of P and gamma, final paragraph
#   Results > Angles (both paragraphs)
#   Figure 3b (the fitted line comes from angle_full)
#
# Numbers this script must produce:
#   ANGLE MODEL
#     intercept  48.5 deg, 95% CI [39.3, 57.8], z = 10.33, P < 0.001
#     slope     -20.2 deg per unit eigenvalue, 95% CI [-34.3, -6.1],
#                z = -2.80, P = 0.005
#     sigma2_study = 203.3, sigma2_obs = 281.2, QE(73) = 5185.7
#     CR2: intercept t(20.1) = 10.01, P < 0.001
#          slope     t(3.02) = -4.67, P = 0.018
#   SIMILARITY MODEL
#     intercept -0.81 on the logit scale, 95% CI [-1.28, -0.35],
#                z = -3.41, P = 0.0007, back-transformed 0.31 [0.22, 0.41]
#     slope      0.83, 95% CI [0.07, 1.60], z = 2.13, P = 0.033
#     sigma2_obs = 0.98, sigma2_study = 0.14, QE(73) = 1745.4
#     CR2: intercept t(15.2) = -3.32, P = 0.005
#          slope     t(2.22) =  1.09, P = 0.38
#
# Requires: df (01), all_angles (13)
# Produces: meta_df, angle_full, sim_full
#
# THE PREDICTOR IS Eigenvalue_1
#   Both reported models use Eigenvalue_1, the raw curvature along the
#   leading Krzanowski axis, not the bias-corrected Eigenvalue_cor_1.
#   The original 08_c_curvature_models.R fitted the similarity model
#   twice, once with each, overwrote the object, and printed both under
#   the same heading. Only the Eigenvalue_1 versions are kept here. If
#   the manuscript is meant to report the bias-corrected predictor, that
#   is a decision to make before rerunning, not a silent swap.
#
# WHY THE ANGLE TRACKS SELECTION MORE SHARPLY THAN THE SIMILARITY
#   Matched directions are ordered by alignment. Averaging over the full
#   subspace dilutes the strong leading-axis signal with the weaker
#   alignment of subordinate directions, so the similarity is expected
#   to respond less sharply than the smallest principal angle, which
#   isolates the best-aligned pair. That expectation is borne out by the
#   CR2 results.
# ======================================================================


# -- 1. Study-level moderators -----------------------------------------

moderators <- df %>%
  dplyr::filter(P.matrix.available. == "Yes") %>%
  dplyr::mutate(Grouping = ifelse(Grouping == "N/A" | Grouping == "",
                                  NA_character_, Grouping)) %>%
  dplyr::distinct(ID, Grouping, Fitness.measurement, .keep_all = TRUE) %>%
  dplyr::transmute(
    StudyID           = as.numeric(ID),
    Grouping,
    Fitness           = as.character(Fitness.measurement),
    sample_size       = as.numeric(Sample.size),
    Setting           = as.character(Setting.field.lab.artificial),
    Type.of.selection = as.character(Type.of.selection),
    Latin.name        = as.character(Latin.name),
    Kingdom           = as.character(Kingdom),
    Class             = as.character(Class)
  )


# -- 2. Merge onto all_angles ------------------------------------------

meta_df <- all_angles %>%
  dplyr::mutate(
    StudyID  = as.numeric(StudyID),
    Grouping = ifelse(Grouping == "N/A" | Grouping == "",
                      NA_character_, Grouping)) %>%
  dplyr::left_join(moderators, by = c("StudyID", "Grouping", "Fitness")) %>%
  dplyr::group_by(StudyID) %>%
  tidyr::fill(sample_size, Setting, Type.of.selection, Latin.name,
              Kingdom, Class, .direction = "downup") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    obs_id   = seq_len(dplyr::n()),
    k_binary = factor(ifelse(Parameters == 1, "k=1", "k>=2"),
                      levels = c("k=1", "k>=2")),
    Setting           = factor(Setting),
    Type.of.selection = factor(Type.of.selection),
    Selection_type    = factor(Selection_type),
    Kingdom           = factor(Kingdom)
  )

relevel_if <- function(f, ref) if (ref %in% levels(f)) stats::relevel(f, ref) else f
meta_df$Setting           <- relevel_if(meta_df$Setting, "Field")
meta_df$Type.of.selection <- relevel_if(meta_df$Type.of.selection, "Natural selection")
meta_df$Selection_type    <- relevel_if(meta_df$Selection_type, "stabilizing")


# -- 3. Logit-transform the similarity ---------------------------------
# Similarity is bounded in [0, 1]. A Smithson-Verkuilen squeeze keeps
# the logit finite at the boundaries; the delta method carries the
# bootstrap variance onto the logit scale.

lg <- logit_squeeze(meta_df$Similarity, meta_df$Similarity_var,
                    N = nrow(meta_df))
meta_df$Similarity_logit     <- lg$logit
meta_df$Similarity_logit_var <- lg$var

stopifnot(all(is.finite(meta_df$Similarity_logit)),
          all(meta_df$Similarity_logit_var > 0, na.rm = TRUE))


# -- 4. Diagnostics ----------------------------------------------------

cat("\n---- 14_alignment_models ----\n")
cat("Comparisons :", nrow(meta_df), "\n")
cat("Studies     :", dplyr::n_distinct(meta_df$StudyID), "\n")
cat("  k = 1  :", sum(meta_df$Parameters == 1, na.rm = TRUE), "\n")
cat("  k >= 2 :", sum(meta_df$Parameters >= 2, na.rm = TRUE), "\n")

cat("\nMissing values per moderator (rows rma.mv will drop):\n")
print(sapply(meta_df[, c("Setting", "Type.of.selection", "Selection_type",
                         "Kingdom", "Eigenvalue_1", "Eigenvalue_cor_1")],
             function(x) sum(is.na(x))))

# The two responses measure related things, so check how related.
ok <- !is.na(meta_df$Angle) & !is.na(meta_df$Similarity)
cat("\nAngle vs Similarity, Pearson r:",
    round(stats::cor(meta_df$Angle[ok], meta_df$Similarity[ok]), 3), "\n")
for (lev in levels(meta_df$k_binary)) {
  sub <- meta_df[meta_df$k_binary == lev & ok, ]
  if (nrow(sub) > 2)
    cat(sprintf("  %-6s r = %.3f (n = %d)\n", lev,
                stats::cor(sub$Angle, sub$Similarity), nrow(sub)))
}

RE2_ang <- list(~ 1 | StudyID, ~ 1 | obs_id)


# ======================================================================
# 5. ANGLE MODEL (degrees)
# ======================================================================
# The smallest principal angle summarises only the best-aligned pair of
# directions. Its intercept is the expected angle at zero correlational
# selection; 45 degrees is what random orientation would give.

angle_full <- metafor::rma.mv(
  yi = Angle, V = Angle_var, mods = ~ Eigenvalue_1,
  random = RE2_ang, data = meta_df, method = "REML")

cat("\n========== ANGLE MODEL (REPORTED) ==========\n")
print(summary(angle_full))

cat("\nIntercept (angle at zero selection):",
    round(as.numeric(angle_full$beta[1]), 1), "degrees, 95% CI [",
    round(angle_full$ci.lb[1], 1), ",", round(angle_full$ci.ub[1], 1), "]\n")
cat("Slope (degrees per unit eigenvalue):",
    round(as.numeric(angle_full$beta[2]), 1), ", 95% CI [",
    round(angle_full$ci.lb[2], 1), ",", round(angle_full$ci.ub[2], 1), "]\n")

cat("\nVariance components:\n")
variance_partition(angle_full, c("study", "observation"))
cat("QE =", round(angle_full$QE, 1), "on",
    angle_full$k - angle_full$p, "df, P =",
    format.pval(angle_full$QEp, digits = 3), "\n")

# Cluster-robust inference with study as the clustering unit.
# The Satterthwaite df for the slope is about 3, which is the diagnostic
# to state once: the slope rests on few studies. It is not a reason to
# hedge the finding itself, which holds under CR2.
cat("\n--- Angle: CR2 cluster-robust (cluster = StudyID) ---\n")
angle_cr2 <- clubSandwich::coef_test(
  angle_full, vcov = "CR2", cluster = angle_full$mf.r[[1]]$StudyID)
print(angle_cr2)


# ======================================================================
# 6. SIMILARITY MODEL (logit scale)
# ======================================================================
# The overall similarity measures overlap across the whole leading
# subspace rather than its best-aligned pair.

sim_full <- metafor::rma.mv(
  yi = Similarity_logit, V = Similarity_logit_var, mods = ~ Eigenvalue_1,
  random = RE2_ang, data = meta_df, method = "REML")

cat("\n========== SIMILARITY MODEL (REPORTED) ==========\n")
print(summary(sim_full))

inv_logit <- function(x) 1 / (1 + exp(-x))
cat("\nIntercept on the logit scale:", round(as.numeric(sim_full$beta[1]), 2),
    ", 95% CI [", round(sim_full$ci.lb[1], 2), ",",
    round(sim_full$ci.ub[1], 2), "]\n")
cat("Back-transformed similarity :", round(inv_logit(sim_full$beta[1]), 2),
    ", 95% CI [", round(inv_logit(sim_full$ci.lb[1]), 2), ",",
    round(inv_logit(sim_full$ci.ub[1]), 2), "]\n")
cat("Slope on the logit scale    :", round(as.numeric(sim_full$beta[2]), 2),
    ", 95% CI [", round(sim_full$ci.lb[2], 2), ",",
    round(sim_full$ci.ub[2], 2), "]\n")

cat("\nVariance components:\n")
variance_partition(sim_full, c("study", "observation"))
cat("QE =", round(sim_full$QE, 1), "on", sim_full$k - sim_full$p,
    "df, P =", format.pval(sim_full$QEp, digits = 3), "\n")

cat("\n--- Similarity: CR2 cluster-robust (cluster = StudyID) ---\n")
sim_cr2 <- clubSandwich::coef_test(
  sim_full, vcov = "CR2", cluster = sim_full$mf.r[[1]]$StudyID)
print(sim_cr2)


# ======================================================================
# 7. SUMMARY TABLE
# ======================================================================

alignment_table <- data.frame(
  model = c("Angle", "Angle", "Similarity (logit)", "Similarity (logit)"),
  term  = c("intercept", "Eigenvalue_1", "intercept", "Eigenvalue_1"),
  estimate = round(c(angle_full$beta, sim_full$beta), 3),
  ci_lower = round(c(angle_full$ci.lb, sim_full$ci.lb), 3),
  ci_upper = round(c(angle_full$ci.ub, sim_full$ci.ub), 3),
  z        = round(c(angle_full$zval, sim_full$zval), 2),
  p        = signif(c(angle_full$pval, sim_full$pval), 3),
  cr2_t    = round(c(angle_cr2$tstat, sim_cr2$tstat), 2),
  cr2_df   = round(c(angle_cr2$df_Satt, sim_cr2$df_Satt), 2),
  cr2_p    = signif(c(angle_cr2$p_Satt, sim_cr2$p_Satt), 3),
  row.names = NULL
)

cat("\n========== ALIGNMENT SUMMARY ==========\n")
print(alignment_table, row.names = FALSE)

utils::write.csv(alignment_table,
                 file.path(PATHS$tables, "alignment_models.csv"),
                 row.names = FALSE)

save_model(list(angle = angle_full, similarity = sim_full,
                angle_cr2 = angle_cr2, similarity_cr2 = sim_cr2,
                meta_df = meta_df, table = alignment_table),
           "14_alignment_models")
