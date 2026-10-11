# ======================================================================
# 04_lande_arnold_raw.R - RAW-DATA STUDIES: LANDE-ARNOLD REFITS AND M
#
# Manuscript sections served:
#   Methods > Effect size calculation
#     "For four studies in which non-linear gradients were not reported
#      but raw trait and fitness data were available, I fitted the Lande
#      and Arnold non-linear selection regression"
#   Methods > Meta-analytic model
#     "Where raw data were available, and multiple gradients were
#      estimated from the same regression within a study, off-diagonal
#      covariances between sampling errors were computed from the
#      variance-covariance matrix of the original regression, yielding a
#      block-diagonal sampling variance matrix V"
#
# Requires: 00_functions.R (run_lande_arnold, compute_study_M,
#           assemble_V_for_metafor), raw files in data/raw_studies/
# Produces: study*_M objects, full
#
# STRUCTURAL RULE
#   compute_study_M() runs every entry of model_specs on the SAME data
#   argument. For sub-models fitted to DIFFERENT individuals (different
#   treatments, replicates, sites), pre-subset the data and make one
#   compute_study_M() call per subset. Only sub-models sharing the same
#   individuals belong in one call with same_individuals = TRUE.
#
# RUNTIME
#   The same_individuals = TRUE branch bootstraps. At B = N_BOOT (1000)
#   this takes minutes, not hours; the original scripts hard-coded
#   B = 10000, which took hours and is not required by the Methods.
# ======================================================================

raw_dir <- PATHS$raw_studies
rp <- function(f) file.path(raw_dir, f)

if (!dir.exists(raw_dir))
  stop("Raw study files not found in ", raw_dir)


# ======================================================================
# STUDY 4 - Calopteryx splendens
# One frame holds FRONT and HIND wing rows for the same individuals.
# The RW1-3 columns are duplicated into front-only and hind-only
# versions so all three wing models run on the pooled frame with
# same_individuals = TRUE.
# ======================================================================

study4_raw <- utils::read.table(
  rp("2008+Klinga+GMM+data+Dryad.txt"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
study4_raw$fitness <- ifelse(study4_raw$Fitness == "Live", 1, 0)

for (i in 1:3) {
  study4_raw[[paste0("RW", i, "_f")]] <-
    ifelse(study4_raw$Wing == "FRONT", study4_raw[[paste0("RW", i)]], NA)
  study4_raw[[paste0("RW", i, "_h")]] <-
    ifelse(study4_raw$Wing == "HIND",  study4_raw[[paste0("RW", i)]], NA)
}

study4_M <- compute_study_M(
  data = study4_raw,
  model_specs = list(
    list(trait_cols  = c("Wing_Length", "Wing_Width",
                         "Patch_Opacity", "Wing_Patch_Centroid_Size"),
         fitness_col = "fitness", binary = TRUE, label = "Both_wings"),
    list(trait_cols  = c("RW1_f", "RW2_f", "RW3_f", "Wing_Width",
                         "Patch_Opacity", "Wing_Patch_Centroid_Size"),
         fitness_col = "fitness", binary = TRUE, label = "Front_wings"),
    list(trait_cols  = c("RW1_h", "RW2_h", "RW3_h", "Wing_Width",
                         "Patch_Opacity", "Wing_Patch_Centroid_Size"),
         fitness_col = "fitness", binary = TRUE, label = "Hind_wings")
  ),
  same_individuals = TRUE, B = N_BOOT, study_id = 4)


# ======================================================================
# STUDY 6 - Poecilia reticulata, four treatments
# Different individuals per treatment, so one call per subset.
# ======================================================================

study6_raw <- data.frame(readr::read_delim(
  rp("MeanChoiceAndColorData.CSV"),
  delim = ";", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE))

traits6 <- c("BodyArea", "Tail", "Orange", "Black",
             "FuzzyBlack", "Silver", "Green", "Purple")

study6_Ms <- lapply(1:4, function(tt) {
  compute_study_M(
    data = study6_raw[study6_raw$Treatment == tt, ],
    model_specs = list(list(trait_cols = traits6, fitness_col = "MeanPref",
                            binary = FALSE, label = paste0("Treatment_", tt))),
    same_individuals = FALSE, study_id = paste0("6_t", tt))
})
study6_t1_M <- study6_Ms[[1]]; study6_t2_M <- study6_Ms[[2]]
study6_t3_M <- study6_Ms[[3]]; study6_t4_M <- study6_Ms[[4]]


# ======================================================================
# STUDY 21 - Platanthera bifolia, control and hand-pollinated
# Treatments applied to different plants, so same_individuals = FALSE:
# M is block-diagonal with no cross-model covariance.
# ======================================================================

study21_raw <- readxl::read_excel(
  rp("Trunschke et al._EcolEvol_2018_PlantHtSpurLSelection_female fitness.xlsx"),
  sheet = "Trunschke et al._raw data")

study21_M <- compute_study_M(
  data = data.frame(study21_raw),
  model_specs = list(
    list(trait_cols  = c("PlantHt", "NbFlrs", "FlSize", "SpurL"),
         fitness_col = "Fitness", binary = FALSE, label = "Control"),
    list(trait_cols  = c("PlantHt", "NbFlrs", "FlSize", "SpurL"),
         fitness_col = "Fitness", binary = FALSE, label = "Hand_pollinated")
  ),
  same_individuals = FALSE, study_id = 21)


# ======================================================================
# STUDY 23 - ant and bird seed dispersal
# Same plants measured for two fitness outcomes.
# ======================================================================

study23_raw <- data.frame(readr::read_csv(rp("Celtis_piecewiseSEM.csv"),
                                          show_col_types = FALSE))

study23_M <- compute_study_M(
  data = study23_raw,
  model_specs = list(
    list(trait_cols  = c("mean.sugar.concentration", "mean.seed.weight"),
         fitness_col = "ant.seed.removal", binary = FALSE,
         label = "ant.seed.removal"),
    list(trait_cols  = c("fruit.crop", "mean.sugar.concentration",
                         "mean.seed.weight"),
         fitness_col = "bird.fruit.removal", binary = FALSE,
         label = "bird.fruit.removal")
  ),
  same_individuals = TRUE, B = N_BOOT, study_id = 23)


# ======================================================================
# STUDY 25 - Pyura dalbyi. Single model; fitness = alive / total.
# ======================================================================

study25_raw <- data.frame(readr::read_csv(rp("Chirgwin_and_Monro_2020.csv"),
                                          show_col_types = FALSE))
study25_raw$fitness <- study25_raw$alive / study25_raw$total

study25_M <- compute_study_M(
  data = study25_raw,
  model_specs = list(
    list(trait_cols  = c("embryo", "time", "larva"),
         fitness_col = "fitness", binary = FALSE, label = "Main")),
  same_individuals = FALSE, study_id = 25)


# ======================================================================
# STUDY 29 - Bugula neritina, two experimental runs (different colonies)
# ======================================================================

study29_raw <- data.frame(readxl::read_excel(
  rp("Metabolic_rate_and_selection_data.xlsx")))

study29_M <- compute_study_M(
  data = study29_raw[study29_raw$Experimental.run == 1, ],
  model_specs = list(
    list(trait_cols  = c("MI.MRE", "Original.colony.size..zooids."),
         fitness_col = "Fitness", binary = FALSE,
         label = "Experimental.run.1")),
  same_individuals = FALSE, study_id = "29_rep1")

study29_rep2_M <- compute_study_M(
  data = study29_raw[study29_raw$Experimental.run == 2, ],
  model_specs = list(
    list(trait_cols  = c("MI.MRE", "Original.colony.size..zooids."),
         fitness_col = "Fitness", binary = FALSE,
         label = "Experimental.run.2")),
  same_individuals = FALSE, study_id = "29_rep2")


# ======================================================================
# STUDY 31 - Ipomopsis aggregata, two sites (different plants)
# The raw file encodes missing data as ".".
# ======================================================================

study31_raw <- data.frame(readr::read_csv(
  rp("Selection_by_pollinators_and_predators.csv"), show_col_types = FALSE))

for (col in c("width", "color", "nectar", "seeds"))
  study31_raw[[col]] <- as.numeric(ifelse(study31_raw[[col]] == ".", NA,
                                          study31_raw[[col]]))

study31_agg_M <- compute_study_M(
  data = study31_raw[study31_raw$Site == "agg", ],
  model_specs = list(list(trait_cols = c("width", "color", "nectar"),
                          fitness_col = "seeds", binary = FALSE, label = "Agg")),
  same_individuals = FALSE, study_id = "31_agg")

study31_hyb_M <- compute_study_M(
  data = study31_raw[study31_raw$Site == "hyb", ],
  model_specs = list(list(trait_cols = c("width", "color", "nectar"),
                          fitness_col = "seeds", binary = FALSE, label = "Hyb")),
  same_individuals = FALSE, study_id = "31_hyb")


# ======================================================================
# STUDY 33 - Raphanus raphanistrum
# Same plants scored for male and female fitness, under two alternative
# focal morphological traits (Exsertion, Separation). The two designs
# share a source file; within each, the two fitness measures come from
# the same plants.
# ======================================================================

study33_raw <- data.frame(readxl::read_excel(
  rp("Waterman_et_al_2022_combined_data.xlsx")))

study33_exs_M <- compute_study_M(
  data = study33_raw,
  model_specs = list(
    list(trait_cols  = c("PC1Flwr", "Exsertion", "OvuleNum",
                         "PollenNum", "PollenViability"),
         fitness_col = "RelSeedSired", binary = FALSE, label = "Exsertion_male"),
    list(trait_cols  = c("PC1Flwr", "Exsertion", "OvuleNum",
                         "PollenNum", "PollenViability"),
         fitness_col = "RelSeedNum", binary = FALSE, label = "Exsertion_female")
  ),
  same_individuals = TRUE, B = N_BOOT, study_id = "33_exs")

study33_sep_M <- compute_study_M(
  data = study33_raw,
  model_specs = list(
    list(trait_cols  = c("PC1Flwr", "Separation", "OvuleNum",
                         "PollenNum", "PollenViability"),
         fitness_col = "RelSeedSired", binary = FALSE, label = "Separation_male"),
    list(trait_cols  = c("PC1Flwr", "Separation", "OvuleNum",
                         "PollenNum", "PollenViability"),
         fitness_col = "RelSeedNum", binary = FALSE, label = "Separation_female")
  ),
  same_individuals = TRUE, B = N_BOOT, study_id = "33_sep")


# ======================================================================
# STUDY 40 - Dryophytes chrysoscelis. Single model.
# ======================================================================

study40_raw <- data.frame(readxl::read_excel(
  rp("Multivariate+selection+in+Hyla+chrysoscelis+-+Evolution.xlsx")))

study40_M <- compute_study_M(
  data = study40_raw,
  model_specs = list(
    list(trait_cols  = c("AlternativeCD", "AlternativeCR", "AlternativeDF",
                         "AlternativeRA", "AlternativePR"),
         fitness_col = "PhonotaxisScore", binary = FALSE, label = "Main")),
  same_individuals = FALSE, study_id = 40)


# ======================================================================
# STUDY 43 - Rangifer tarandus. Single model.
# ======================================================================

study43_raw <- utils::read.table(rp("natural_selection_data.txt"),
                                 header = TRUE, sep = " ",
                                 stringsAsFactors = FALSE)

study43_M <- compute_study_M(
  data = study43_raw,
  model_specs = list(
    list(trait_cols  = c("birthdate_sd", "birthmass_sd"),
         fitness_col = "lambda_rel", binary = FALSE, label = "Main")),
  same_individuals = FALSE, study_id = 43)


# ======================================================================
# ASSEMBLE
# Block-diagonal across studies. Studies are independent at the sampling
# level. `full$V` covers only the raw-data studies; 05_build_combined.R
# pads it with a diagonal block for the published-value studies.
# ======================================================================

raw_M_list <- list(
  study4_M,
  study6_t1_M, study6_t2_M, study6_t3_M, study6_t4_M,
  study21_M, study23_M, study25_M,
  study29_M, study29_rep2_M,
  study31_agg_M, study31_hyb_M,
  study33_exs_M, study33_sep_M,
  study40_M, study43_M
)

full <- assemble_V_for_metafor(raw_M_list)

cat("\n---- 04_lande_arnold_raw ----\n")
cat("Sub-models fitted   :", length(raw_M_list), "\n")
cat("Raw-data V dimension:", paste(dim(full$V), collapse = " x "), "\n")
cat("Non-linear gradients:", length(full$gamma_hat), "\n")

save_model(raw_M_list, "raw_study_M_objects")
