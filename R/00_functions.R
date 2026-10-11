# ======================================================================
# 00_functions.R - SHARED HELPER FUNCTIONS
#
# Every function used by more than one analysis script lives here.
#
# ======================================================================


# ======================================================================
# SECTION A - MODEL INFERENCE HELPERS
# ======================================================================

#' Cluster-robust (CR2) summary of an rma.mv fit
#'
#' Prints coefficient-level CR2 tests with Satterthwaite degrees of
#' freedom, robust confidence intervals, and a joint Wald test per
#' multi-level moderator. Used for every inference reported in the
#' Results (Methods: "Inference on predictors used cluster-robust (CR2)
#' standard errors with study as the clustering unit").
rve_summary <- function(model, cluster_vec, label = "") {
  cat("\n========== RVE (CR2) :", label, "==========\n")
  
  # rma.mv silently drops rows with a missing moderator, so a cluster
  # vector taken from the full data frame can be longer than the fitted
  # model. Fall back to the study column in the model's own random-effects
  # frame when the lengths disagree, rather than letting clubSandwich fail
  # on a length mismatch.
  if (length(cluster_vec) != model$k) {
    cluster_vec <- tryCatch(model$mf.r[[1]][[1]],
                            error = function(e) NULL)
    if (is.null(cluster_vec) || length(cluster_vec) != model$k) {
      cat("Cluster vector does not match the fitted model (k = ", model$k,
          "). Skipping CR2.\n", sep = "")
      return(invisible(NULL))
    }
    cat("Note: cluster vector rebuilt from the model frame",
        "(rows were dropped for missing moderators).\n")
  }
  
  ct <- clubSandwich::coef_test(model, vcov = "CR2", cluster = cluster_vec)
  print(ct)
  
  ci <- clubSandwich::conf_int(model, vcov = "CR2", cluster = cluster_vec)
  cat("\n-- Robust 95% CIs --\n"); print(ci)
  
  has_mods <- !is.null(model$formula.mods) &&
    length(attr(stats::terms(model$formula.mods), "term.labels")) > 0
  
  if (has_mods) {
    mod_names  <- attr(stats::terms(model$formula.mods), "term.labels")
    mod_coefs  <- stats::coef(model)
    beta_names <- names(mod_coefs)
    
    cat("\n-- Joint Wald tests per moderator --\n")
    for (m in mod_names) {
      idx <- which(startsWith(beta_names, m))
      if (length(idx) == 0) { cat("\n", m, ": no matching coefficients\n"); next }
      wt <- tryCatch(
        clubSandwich::Wald_test(
          model,
          constraints = clubSandwich::constrain_zero(idx, mod_coefs),
          vcov = "CR2", cluster = cluster_vec),
        error = function(e) paste("skipped:", conditionMessage(e))
      )
      cat("\n", m, ":\n"); print(wt)
    }
  }
  invisible(list(coef = ct, ci = ci))
}


#' Variance partition across the random levels of an rma.mv fit
#'
#' Reports each sigma^2 as a percentage of the total random-effects
#' variance. Results: "42%, 37%, and 20% of the random-effects variance".
variance_partition <- function(model, labels = NULL) {
  s2 <- model$sigma2
  if (is.null(labels)) labels <- paste0("level_", seq_along(s2))
  out <- data.frame(level = labels, sigma2 = s2,
                    pct = 100 * s2 / sum(s2))
  print(out, row.names = FALSE, digits = 4)
  invisible(out)
}


#' Folded-normal expected absolute gradient (Morrissey 2016)
#'
#' E[|gamma|] = tau * sqrt(2/pi) * exp(-mu^2 / 2 tau^2)
#'              + mu * (2 * Phi(mu / tau) - 1)
#' with mu the signed meta-analytic mean and tau^2 the summed
#' random-effect variances of the SIGNED model.
fold_norm <- function(fit) {
  mu     <- as.numeric(fit$beta[1])
  sigma2 <- sum(fit$sigma2)
  sigma  <- sqrt(sigma2)
  sigma * sqrt(2 / pi) * exp(-mu^2 / (2 * sigma2)) +
    mu * (2 * stats::pnorm(mu / sigma) - 1)
}


#' Non-parametric cluster bootstrap of the folded-normal E[|gamma|]
#'
#' Resamples STUDIES with replacement, refits the signed model, and
#' recomputes E[|gamma|]. Percentile intervals are taken across
#' replicates (Methods: B = 1,000).
boot_fold <- function(data, B = N_BOOT, formula_random, seed = SEED) {
  set.seed(seed)
  studies <- unique(data$study_id)
  out     <- numeric(B)
  
  for (b in seq_len(B)) {
    pick <- sample(studies, length(studies), replace = TRUE)
    boot_df <- do.call(rbind, lapply(seq_along(pick), function(i) {
      sub <- data[data$study_id == pick[i], ]
      sub$study_id <- paste0(pick[i], "_rep", i)
      sub$obs_id   <- paste0(sub$obs_id, "_rep", i)
      if ("group_id" %in% names(sub))
        sub$group_id <- paste0(sub$group_id, "_rep", i)
      sub
    }))
    
    fit_b <- tryCatch(
      metafor::rma.mv(yi, vi, random = formula_random,
                      data = boot_df, method = "REML"),
      error = function(e) NULL
    )
    out[b] <- if (is.null(fit_b)) NA_real_ else fold_norm(fit_b)
  }
  out
}


# ======================================================================
# SECTION B - KRZANOWSKI SUBSPACE COMPARISON
# ======================================================================

#' Krzanowski (1979) subspace alignment between gamma and P
#'
#' Builds the gamma and P matrices for one (study x grouping x fitness)
#' subset, compares their leading k = max(1, floor(n/2)) subspaces via
#' S = A'BB'A, and returns the smallest principal angle, the overall
#' similarity sum(lambda)/k, the leading-eigenvector angle, per-axis
#' curvatures, Reynolds-style null p-values and bias-corrected
#' curvatures, and parametric bootstrap SEs.
#'
#' @param subset_df rows of df for one comparison
#' @param n_sim bootstrap and null replicates (Methods: B = 1000)
#'   note in 00_setup.R and README "Known discrepancies".
calc_angle <- function(subset_df, n_sim = N_BOOT) {
  tryCatch({
    
    # Normalise "N/A" strings to real NA
    subset_df$Trait.1[subset_df$Trait.1 == "N/A"] <- NA
    subset_df$Trait.2[subset_df$Trait.2 == "N/A"] <- NA
    
    for (v in c("Gamma.ii", "Gamma.ij", "Gamma.ii.se", "Gamma.ij.se",
                "Pii", "Pij"))
      subset_df[[v]] <- as.numeric(subset_df[[v]])
    
    traits <- unique(c(subset_df$Trait.1, subset_df$Trait.2))
    traits <- traits[!is.na(traits)]
    n <- length(traits)
    if (n < 2) return(NULL)
    
    Gamma <- matrix(0, n, n, dimnames = list(traits, traits))
    P     <- matrix(0, n, n, dimnames = list(traits, traits))
    SE_G  <- matrix(0, n, n, dimnames = list(traits, traits))
    
    diag_rows <- subset_df[is.na(subset_df$Trait.2), ]
    off_rows  <- subset_df[!is.na(subset_df$Trait.2), ]
    
    if (nrow(diag_rows) > 0) {
      i <- match(diag_rows$Trait.1, traits); v <- !is.na(i)
      gii    <- diag_rows$Gamma.ii[v]
      gii_se <- diag_rows$Gamma.ii.se[v]
      Gamma[cbind(i[v], i[v])] <- gii
      P[cbind(i[v], i[v])]     <- diag_rows$Pii[v]
      SE_G[cbind(i[v], i[v])]  <- gii_se
    }
    if (nrow(off_rows) > 0) {
      i <- match(off_rows$Trait.1, traits)
      j <- match(off_rows$Trait.2, traits)
      v <- !is.na(i) & !is.na(j)
      Gamma[cbind(i[v], j[v])] <- Gamma[cbind(j[v], i[v])] <- off_rows$Gamma.ij[v]
      P[cbind(i[v], j[v])]     <- P[cbind(j[v], i[v])]     <- off_rows$Pij[v]
      SE_G[cbind(i[v], j[v])]  <- SE_G[cbind(j[v], i[v])]  <- off_rows$Gamma.ij.se[v]
    }
    SE_G[is.na(SE_G) | is.infinite(SE_G) | is.nan(SE_G)] <- 0
    
    if (any(is.na(Gamma)) || any(is.infinite(Gamma))) return(NULL)
    if (any(is.na(P))     || any(is.infinite(P)))     return(NULL)
    
    Gamma <- (Gamma + t(Gamma)) / 2
    P     <- (P + t(P)) / 2
    
    # Null distribution centred on zero (Reynolds-style)
    null_eigs <- replicate(n_sim, {
      gn <- matrix(stats::rnorm(n * n, mean = 0, sd = SE_G), n, n)
      gn <- (gn + t(gn)) / 2
      abs(eigen(gn, symmetric = TRUE)$values)
    })
    
    eigGamma <- eigen(Gamma, symmetric = TRUE)
    eigP     <- eigen(P,     symmetric = TRUE)
    raw_gamma_eigs <- eigGamma$values
    
    # Leading-eigenvector angle (theta_max)
    topk <- function(eig, k) order(abs(eig$values), decreasing = TRUE)[1:k]
    Gmax <- eigGamma$vectors[, topk(eigGamma, 1)]
    Pmax <- eigP$vectors[, 1]
    angle_max <- acos(sum(Gmax * Pmax) /
                        (sqrt(sum(Gmax^2)) * sqrt(sum(Pmax^2)))) * 180 / pi
    if (angle_max > 90) angle_max <- 180 - angle_max
    
    k <- max(1, floor(n / 2))
    normalize_vecs <- function(m) apply(m, 2, function(v) v / sqrt(sum(v^2)))
    A <- normalize_vecs(eigP$vectors[,     topk(eigP,     k), drop = FALSE])
    B <- normalize_vecs(eigGamma$vectors[, topk(eigGamma, k), drop = FALSE])
    
    S         <- t(A) %*% B %*% t(B) %*% A
    S_decomp  <- eigen(S, symmetric = TRUE)
    S_values  <- S_decomp$values
    V_aligned <- A %*% S_decomp$vectors
    
    # Curvature along each Krzanowski axis
    abs_curvatures <- sapply(seq_len(ncol(V_aligned)), function(i) {
      v <- V_aligned[, i]
      sum(abs(eigGamma$values) * as.vector(t(eigGamma$vectors) %*% v)^2)
    })
    
    eigen_results <- lapply(seq_len(ncol(V_aligned)), function(i) {
      v   <- V_aligned[, i]
      obs <- abs_curvatures[i]
      null_curv <- apply(null_eigs, 2, function(eigs)
        sum(eigs * as.vector(t(eigGamma$vectors) %*% v)^2))
      list(pval = mean(null_curv >= obs), corrected = obs - mean(null_curv))
    })
    eigen_pvals     <- sapply(eigen_results, `[[`, "pval")
    eigen_corrected <- sapply(eigen_results, `[[`, "corrected")
    
    # Selection regime per axis and overall
    curvatures <- apply(V_aligned, 2, function(v) as.numeric(t(v) %*% Gamma %*% v))
    sel_label <- function(x) ifelse(x < 0, "stabilizing",
                                    ifelse(x > 0, "disruptive", "neutral"))
    selection_per_axis <- sel_label(curvatures)
    selection_type1    <- sel_label(curvatures[1])
    selection_overall  <- if (all(curvatures < 0)) "stabilizing"
    else if (all(curvatures > 0)) "disruptive"
    else if (all(curvatures == 0)) "neutral"
    else "saddle point"
    
    S_clip        <- pmax(pmin(S_values, 1), 0)
    angles_sorted <- sort(acos(sqrt(S_clip)) * 180 / pi)
    similarity    <- sum(S_values) / k
    
    # Parametric bootstrap of gamma: SE for BOTH summary statistics.
    # The similarity SE is new - see README "Changes made". Without it
    # the similarity model in 14_alignment_models.R has no sampling
    # variance to weight by.
    boot <- replicate(n_sim, {
      gb <- matrix(stats::rnorm(n * n, mean = as.vector(Gamma),
                                sd = as.vector(SE_G)), n, n)
      gb <- (gb + t(gb)) / 2
      eB <- eigen(gb, symmetric = TRUE)
      Bb <- normalize_vecs(eB$vectors[, topk(eB, k), drop = FALSE])
      Sb <- t(A) %*% Bb %*% t(Bb) %*% A
      sv <- eigen(Sb, symmetric = TRUE)$values
      sb <- pmax(pmin(sv, 1), 0)
      ang <- acos(sqrt(sb)) * 180 / pi
      ang[ang > 90] <- 180 - ang[ang > 90]
      c(min(ang), sum(sv) / k)
    })
    Angle_SE      <- stats::sd(boot[1, ], na.rm = TRUE)
    Similarity_SE <- stats::sd(boot[2, ], na.rm = TRUE)
    
    pad4 <- function(x, prefix)
      stats::setNames(lapply(1:4, function(i) if (i <= length(x)) x[i] else NA),
                      paste0(prefix, 1:4))
    
    out <- tibble::tibble(
      Angle           = if (length(angles_sorted) > 0) min(angles_sorted) else NA_real_,
      Angle_SE        = Angle_SE,
      Angle_var       = Angle_SE^2,
      Angle_max       = angle_max,
      Similarity      = similarity,
      Similarity_SE   = Similarity_SE,
      Similarity_var  = Similarity_SE^2,
      n_traits        = n,
      Parameters      = k,
      Selection_type  = selection_overall,
      Selection_type1 = selection_type1
    )
    cbind(
      out,
      as.data.frame(pad4(angles_sorted,      "Angle_"),           check.names = FALSE),
      as.data.frame(pad4(selection_per_axis, "Axis_"),            check.names = FALSE),
      as.data.frame(pad4(abs_curvatures,     "Eigenvalue_"),      check.names = FALSE),
      as.data.frame(pad4(eigen_pvals,        "Eigenvalue_pval_"), check.names = FALSE),
      as.data.frame(pad4(eigen_corrected,    "Eigenvalue_cor_"),  check.names = FALSE),
      as.data.frame(pad4(raw_gamma_eigs,     "Raw_Gamma_eig_"),   check.names = FALSE)
    )
  }, error = function(e) {
    message("Skipping subset: ", e$message); NULL
  })
}
# ======================================================================
# SECTION C - SMALL UTILITIES
# ======================================================================

#' Smithson-Verkuilen squeeze then logit, with a delta-method variance
#'
#' Similarity is bounded in [0, 1]; squeezing avoids infinite logits at
#' the boundaries. var(logit(s)) = var(s) / (s (1 - s))^2.
logit_squeeze <- function(s, var_s = NULL, N = length(s)) {
  ss <- (s * (N - 1) + 0.5) / N
  out <- list(logit = log(ss / (1 - ss)))
  if (!is.null(var_s)) out$var <- var_s / (ss * (1 - ss))^2
  out
}

#' Significance stars for figure annotation
p_stars <- function(p)
  ifelse(p < .001, "***", ifelse(p < .01, "**",
                                 ifelse(p < .05, "*", ifelse(p < .1, ".", ""))))

# Type-of-selection palette, after Jarrett et al. 2025 (Fig. 2):
# amber yellow and crimson red
pal <- c("Natural selection" = "#F0B429",
         "Sexual selection"  = "#B5172B")

# Text sizes for figures, set for reading from the back of a room
FIG_BASE_PT  <- 20   # theme base size, in points
FIG_ANNOT_SZ <- 6    # geom_text size, in mm (about 17 pt)

#' Save a fitted model object to output/models/
save_model <- function(obj, name) {
  saveRDS(obj, file.path(PATHS$models, paste0(name, ".rds")))
  invisible(obj)
}