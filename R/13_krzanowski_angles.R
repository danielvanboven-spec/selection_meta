# ======================================================================
# 13_krzanowski_angles.R - SUBSPACE COMPARISON OF P AND gamma
#
# Manuscript sections served:
#   Methods > Subspace comparison of P and gamma
#   Results > Angles (the 75 comparisons from 26 studies)
#   Figure 3
#
# Requires: df (01_load_clean.R), calc_angle() (00_functions.R)
# Produces: all_angles, output/tables/all_angles.csv
#
# WHAT calc_angle() RETURNS PER COMPARISON
#   Angle           smallest Krzanowski principal angle, theta_min
#   Angle_SE/_var   parametric bootstrap SE and variance of theta_min
#   Angle_max       angle between the leading eigenvectors of gamma and P
#   Similarity      overall subspace similarity, sum(lambda_i) / k
#   Similarity_SE   bootstrap SE of the similarity (NEW, see below)
#   Parameters      subspace dimension k = max(1, floor(n / 2))
#   n_traits        number of traits n in the comparison
#   Eigenvalue_*    curvature along each Krzanowski axis
#   Eigenvalue_cor_* the same, bias-corrected against the null
#   Selection_type  stabilising, disruptive, saddle point or neutral
#
# TWO FIXES RELATIVE TO THE ORIGINAL 08_meta_curvature.R
#   1. The original had a stray "W" after compute_angles_for_study(df,
#      sid), which is a syntax error: the script could not be sourced.
#   2. The bootstrap recorded only Angle_SE. 14_alignment_models.R needs
#      a sampling variance for the similarity as well, and referred to a
#      column Similarity_logit_var that nothing ever created. calc_angle()
#      now returns Similarity_SE from the same bootstrap replicates, and
#      14 derives the logit variance from it by the delta method.
#
# k = max(1, floor(n / 2))
#   Krzanowski's construction requires k <= n / 2. For n <= 3 that
#   forces k = 1, which reduces the comparison to the single leading
#   eigenvector and returns Angle == Angle_max. Only comparisons with
#   n >= 4 traits give a genuine subspace angle. There are 24 of those,
#   which is the limitation discussed in the manuscript.
#
# RUNTIME
#   Two sets of N_BOOT replicates per comparison. Expect a few minutes.
# ======================================================================


# -- 1. Loop over study x grouping x fitness measure -------------------
# Only studies that report a phenotypic covariance matrix can enter.
# Rows whose Grouping is "N/A" belong to every grouping within a study,
# so they are pooled in with each one.

compute_angles_for_study <- function(df, studyID) {

  study <- df[df$ID == studyID & df$P.matrix.available. == "Yes", ]
  if (nrow(study) == 0) return(NULL)

  groups    <- unique(study$Grouping)
  fitnesses <- unique(study$Fitness.measurement)
  results   <- list()

  for (g in groups) {
    if (is.na(g) || g == "N/A" || g == "") next

    fitnesses_g <- unique(study$Fitness.measurement[
      study$Grouping == g | study$Grouping == "N/A"])

    for (f in fitnesses_g) {
      if (is.na(f) || f == "N/A" || f == "") next

      has_group <- any(study$Grouping == g & study$Fitness.measurement == f,
                       na.rm = TRUE)
      if (!has_group) next

      subset_df <- study[(study$Grouping == g | study$Grouping == "N/A") &
                           study$Fitness.measurement == f, ]
      if (nrow(subset_df) == 0) next

      res <- calc_angle(subset_df)
      if (is.null(res) || nrow(res) == 0) next

      results[[length(results) + 1]] <- cbind(
        data.frame(StudyID = studyID, Grouping = g, Fitness = f,
                   stringsAsFactors = FALSE), res)
    }
  }

  # Fallback for studies with no usable grouping: loop fitness only
  if (length(results) == 0) {
    for (f in fitnesses) {
      if (is.na(f) || f == "N/A" || f == "") next
      subset_df <- study[study$Fitness.measurement == f, ]
      if (nrow(subset_df) == 0) next

      res <- calc_angle(subset_df)
      if (is.null(res) || nrow(res) == 0) next

      results[[length(results) + 1]] <- cbind(
        data.frame(StudyID = studyID, Grouping = NA_character_, Fitness = f,
                   stringsAsFactors = FALSE), res)
    }
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}


# -- 2. Run ------------------------------------------------------------

set.seed(SEED)

study_ids <- sort(stats::na.omit(unique(
  df$ID[df$P.matrix.available. == "Yes" & !is.na(df$ID)])))

cat("\n---- 13_krzanowski_angles ----\n")
cat("Studies reporting a P matrix:", length(study_ids), "\n")

results_list <- lapply(study_ids, function(sid) {
  cat("  study", sid, "\n")
  compute_angles_for_study(df, sid)
})

all_angles <- do.call(rbind, results_list)


# -- 3. Report ---------------------------------------------------------

cat("\n---- Comparisons ----\n")
cat("Total gamma-P comparisons :", nrow(all_angles), " (Results: 75)\n")
cat("Contributing studies      :", dplyr::n_distinct(all_angles$StudyID),
    " (Results: 26)\n")
cat("  k = 1 :", sum(all_angles$Parameters == 1), "\n")
cat("  k = 2 :", sum(all_angles$Parameters == 2), " (Results: 24)\n")
cat("  k >= 3:", sum(all_angles$Parameters >= 3), "\n")

cat("\nSmallest principal angle : mean",
    round(mean(all_angles$Angle, na.rm = TRUE), 1), "degrees\n")
cat("Leading-eigenvector angle: mean",
    round(mean(all_angles$Angle_max, na.rm = TRUE), 1), "degrees\n")

# Gain from the subspace comparison over the single-vector one.
# Figure 3a: "the mean gain across k >= 2 comparisons is 45.7 degrees".
# The gain is zero by construction at k = 1, so only k >= 2 counts.
delta_k2 <- with(all_angles[all_angles$Parameters >= 2, ], Angle_max - Angle)
cat("Mean gain at k >= 2      :", round(mean(delta_k2, na.rm = TRUE), 1),
    "degrees over", length(delta_k2), "comparisons\n")

cat("\nSelection regime:\n")
print(table(all_angles$Selection_type, useNA = "ifany"))

utils::write.csv(all_angles,
                 file.path(PATHS$tables, "all_angles.csv"), row.names = FALSE)
save_model(all_angles, "13_all_angles")

cat("\nWritten to output/tables/all_angles.csv\n")
