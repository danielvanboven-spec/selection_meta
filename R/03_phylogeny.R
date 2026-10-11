# ======================================================================
# 03_phylogeny.R - PHYLOGENETIC TREE
#
# Manuscript sections served:
#   Methods > Statistical analysis > Phylogenetic tree
#   Figure 1
#
# Requires: df (01_load_clean.R)
# Produces: tree, tree_bl, phylo_cor_matrix, tip_classes, phylo_standins
#           output/tables/phylogeny_standins.csv
#
# The decorated version of the tree (class colours, phylopics, photo
# panel) is figures/fig01_phylogeny.R. This script only builds the tree
# object and its correlation matrix.
# ======================================================================

if (!load_optional(pkgs_phylo))
  stop("03_phylogeny.R needs rotl, ape and ggtree.")

# Optional manual stand-ins, used instead of the automatic search.
# Format: "name in df" = "species to place it by". Leave empty to let
# the search decide.
MANUAL_STANDIN <- c()

# How far up the taxonomy the search may go when there is no family
# among the ancestors
MAX_LEVELS_UP <- 6


# -- 1. Species list ---------------------------------------------------

TreeSpecies <- unique(df[, c("Kingdom", "Class", "Latin.name")])
TreeSpecies <- TreeSpecies[!is.na(TreeSpecies$Latin.name) &
                             trimws(TreeSpecies$Latin.name) != "" &
                             TreeSpecies$Latin.name != "N/A", ]

dup <- unique(TreeSpecies$Latin.name[duplicated(TreeSpecies$Latin.name)])
if (length(dup) > 0) {
  warning("Species with more than one Kingdom/Class in df (first kept): ",
          paste(dup, collapse = ", "), call. = FALSE)
  TreeSpecies <- TreeSpecies[!duplicated(TreeSpecies$Latin.name), ]
}

sp <- TreeSpecies
rownames(sp) <- NULL

cat("\n---- Phylogeny ----\n")
cat("Distinct species in df:", nrow(sp), "\n")


# -- 2. Match to the Open Tree of Life ---------------------------------
# "Genus spp." is recorded to genus only, so the genus is what is looked
# up and a species inside it is used as the stand-in.

sp$genus_level <- grepl("\\s+spp?\\.?$", sp$Latin.name)
sp$query       <- trimws(sub("\\s+spp?\\.?$", "", sp$Latin.name))

m <- rotl::tnrs_match_names(sp$query)
m <- m[match(tolower(sp$query), tolower(m$search_string)), ]

sp$ott_id   <- as.numeric(m$ott_id)
sp$otl_name <- m$unique_name

# The Open Tree taxonomy contains wrong synonyms: it returns the yam
# Dioscorea laxiflora for the hackberry Celtis ehrenbergiana. A real
# synonym keeps the genus or the species epithet (Hyla chrysoscelis for
# Dryophytes chrysoscelis), so a match that changes both is rejected and
# the species is placed by a stand-in from its own genus instead.
word_of <- function(x, k) vapply(strsplit(trimws(as.character(x)), "\\s+"),
                                 function(p) if (length(p) >= k) p[k] else NA_character_,
                                 character(1))
stem_of <- function(x) substr(tolower(x), 1, 5)

same_genus   <- !is.na(sp$otl_name) &
  tolower(word_of(sp$query, 1)) == tolower(word_of(sp$otl_name, 1))
same_epithet <- !sp$genus_level &
  !is.na(word_of(sp$query, 2)) & !is.na(word_of(sp$otl_name, 2)) &
  stem_of(word_of(sp$query, 2)) == stem_of(word_of(sp$otl_name, 2))

rejected <- !is.na(sp$ott_id) & !same_genus & !same_epithet
if (any(rejected)) {
  cat("Open Tree matches rejected (neither genus nor epithet agrees):\n")
  print(data.frame(species_in_df = sp$Latin.name[rejected],
                   open_tree_match = sp$otl_name[rejected]),
        row.names = FALSE)
  sp$ott_id[rejected]   <- NA_real_
  sp$otl_name[rejected] <- NA_character_
}

has_id <- !is.na(sp$ott_id)
sp$in_tree <- FALSE
sp$in_tree[has_id] <- unname(rotl::is_in_tree(sp$ott_id[has_id]))

# Species the synthetic tree places directly
direct <- has_id & sp$in_tree & !sp$genus_level

sp$tree_ott  <- ifelse(direct, sp$ott_id, NA_real_)
sp$tree_name <- ifelse(direct, sp$otl_name, NA_character_)
sp$placed_by <- ifelse(direct, "species itself", NA_character_)

# Two names in df that resolve to one Open Tree taxon are one species
clash <- sp$tree_ott[direct][duplicated(sp$tree_ott[direct])]
if (length(clash) > 0)
  stop("Different names in df resolve to the same Open Tree taxon: ",
       paste(sp$Latin.name[sp$tree_ott %in% clash], collapse = ", "),
       "\n  Harmonise them in 01_load_clean.R.", call. = FALSE)


# -- 3. Stand-ins for species the synthetic tree cannot place ----------

# A taxon and its ancestors, nearest first
otl_nodes <- function(ott_id, include_self) {
  info <- rotl::taxonomy_taxon_info(ott_id, include_lineage = TRUE)[[1]]
  lin  <- info$lineage
  out <- data.frame(
    ott_id = vapply(lin, function(z) as.numeric(z$ott_id), numeric(1)),
    name   = vapply(lin, function(z) as.character(z$name), character(1)),
    rank   = vapply(lin, function(z) as.character(z$rank), character(1)),
    stringsAsFactors = FALSE
  )
  if (include_self)
    out <- rbind(data.frame(ott_id = as.numeric(ott_id),
                            name   = as.character(info$name),
                            rank   = as.character(info$rank),
                            stringsAsFactors = FALSE), out)
  
  # Do not search above the family
  fam <- which(out$rank == "family")
  if (length(fam) > 0) out[seq_len(fam[1]), ] else utils::head(out, MAX_LEVELS_UP)
}

# Labels of the form Genus_species_ott123 to a table of plain binomials
# (no "sp.", no hybrids, no subspecies), in alphabetical order
labels_to_species <- function(lab) {
  lab <- as.character(unlist(lab))
  lab <- lab[grepl("_ott[0-9]+$", lab)]
  out <- data.frame(
    name   = gsub("_", " ", sub("_ott[0-9]+$", "", lab)),
    ott_id = as.numeric(sub("^.*_ott([0-9]+)$", "\\1", lab)),
    stringsAsFactors = FALSE
  )
  out <- out[grepl("^[A-Z][a-z]+ [a-z][a-z-]+$", out$name), ]
  out <- out[!duplicated(out$ott_id), ]
  out[order(out$name), ]
}

# First species below a taxon that is in the synthetic tree and not used.
#
# Route 1: the taxon is a node of the synthetic tree, so its subtree lists
#   the species directly.
# Route 2: the taxon is "broken" (not one clade in the synthetic tree), so
#   it has no subtree. Its species are then taken from the taxonomy and
#   checked one by one. This is the case for Celtis.
MAX_CHECKS <- 300   # most species checked one by one under route 2

first_species_in_tree <- function(ott_id, used) {
  
  sub_tree <- tryCatch(rotl::tol_subtree(ott_id = ott_id),
                       error = function(e) NULL)
  if (!is.null(sub_tree)) {
    tips <- labels_to_species(sub_tree$tip.label)
    tips <- tips[!tips$ott_id %in% used, ]
    if (nrow(tips) > 0) return(tips[1, ])
  }
  
  taxa <- tryCatch(rotl::taxonomy_subtree(ott_id, output_format = "taxa"),
                   error = function(e) NULL)
  if (is.null(taxa)) return(NULL)
  cand <- labels_to_species(taxa)
  cand <- cand[!cand$ott_id %in% used, ]
  for (j in seq_len(min(nrow(cand), MAX_CHECKS))) {
    ok <- tryCatch(isTRUE(unname(rotl::is_in_tree(cand$ott_id[j]))),
                   error = function(e) FALSE)
    if (ok) return(cand[j, ])
  }
  NULL
}

find_standin <- function(start_id, include_self, used) {
  nodes <- otl_nodes(start_id, include_self)
  for (i in seq_len(nrow(nodes))) {
    cat("    searching", nodes$rank[i], nodes$name[i], "... ")
    hit <- first_species_in_tree(nodes$ott_id[i], used)
    if (is.null(hit)) { cat("nothing usable\n"); next }
    cat("found", hit$name, "\n")
    return(list(ott_id = hit$ott_id, name = hit$name,
                via = paste(nodes$rank[i], nodes$name[i])))
  }
  NULL
}

for (i in which(is.na(sp$tree_ott))) {
  
  used <- sp$tree_ott[!is.na(sp$tree_ott)]
  nm   <- sp$Latin.name[i]
  cat("  Stand-in for", nm, "\n")
  
  # Manual stand-in, if one was given
  if (nm %in% names(MANUAL_STANDIN)) {
    mm <- rotl::tnrs_match_names(MANUAL_STANDIN[[nm]])
    ok <- !is.na(mm$ott_id[1]) && isTRUE(unname(rotl::is_in_tree(mm$ott_id[1])))
    if (!ok)
      stop("Manual stand-in for ", nm, " is not in the synthetic tree: ",
           MANUAL_STANDIN[[nm]], call. = FALSE)
    if (mm$ott_id[1] %in% used)
      stop("Manual stand-in for ", nm, " is already a tip: ",
           MANUAL_STANDIN[[nm]], call. = FALSE)
    sp$tree_ott[i]  <- as.numeric(mm$ott_id[1])
    sp$tree_name[i] <- mm$unique_name[1]
    sp$placed_by[i] <- "manual stand-in"
    next
  }
  
  if (!is.na(sp$ott_id[i])) {
    # Known to the taxonomy. A genus is searched from the genus itself,
    # a species from its parent upwards.
    hit <- find_standin(sp$ott_id[i], include_self = sp$genus_level[i], used)
  } else {
    # Unknown to the taxonomy, or its match was rejected: fall back on
    # the genus name
    genus <- strsplit(sp$query[i], "\\s+")[[1]][1]
    g <- tryCatch(rotl::tnrs_match_names(genus), error = function(e) NULL)
    genus_ok <- !is.null(g) && !is.na(g$ott_id[1]) &&
      tolower(word_of(g$unique_name[1], 1)) == tolower(genus)
    if (!genus_ok)
      cat("    genus", genus, "not found in the Open Tree taxonomy",
          if (!is.null(g)) paste0("(matched: ", g$unique_name[1], ")"), "\n")
    hit <- if (!genus_ok) NULL else
      find_standin(as.numeric(g$ott_id[1]), include_self = TRUE, used)
  }
  
  if (is.null(hit)) next
  
  sp$tree_ott[i]  <- hit$ott_id
  sp$tree_name[i] <- hit$name
  sp$placed_by[i] <- paste("stand-in from", hit$via)
}

unplaced <- sp$Latin.name[is.na(sp$tree_ott)]
if (length(unplaced) > 0)
  stop("No stand-in found for: ", paste(unplaced, collapse = ", "),
       "\n  Add one to MANUAL_STANDIN at the top of this script.",
       call. = FALSE)

phylo_standins <- sp[sp$placed_by != "species itself",
                     c("Latin.name", "tree_name", "placed_by")]
names(phylo_standins) <- c("species_in_df", "placed_by_species", "found_at")
rownames(phylo_standins) <- NULL

cat("Placed directly        :", sum(sp$placed_by == "species itself"), "\n")
cat("Placed by a stand-in   :", nrow(phylo_standins), "\n")
if (nrow(phylo_standins) > 0) print(phylo_standins, row.names = FALSE)

utils::write.csv(phylo_standins,
                 file.path(PATHS$tables, "phylogeny_standins.csv"),
                 row.names = FALSE)


# -- 4. Induced subtree ------------------------------------------------
# Tips are matched back to df by Open Tree id, not by name, so synonyms
# and stand-ins both end up labelled with the name used in df.

tree <- rotl::tol_induced_subtree(ott_ids = sp$tree_ott)

tip_ott <- suppressWarnings(
  as.numeric(sub("^.*ott([0-9]+)$", "\\1", tree$tip.label)))
idx <- match(tip_ott, sp$tree_ott)

if (anyNA(idx))
  stop("The Open Tree returned tips that match no requested taxon: ",
       paste(tree$tip.label[is.na(idx)], collapse = ", "), call. = FALSE)

tree$tip.label <- sp$Latin.name[idx]

not_in_tree <- setdiff(sp$Latin.name, tree$tip.label)
if (length(not_in_tree) > 0)
  stop("Requested but missing from the returned tree: ",
       paste(not_in_tree, collapse = ", "), call. = FALSE)

cat("Tree tips:", length(tree$tip.label), "\n")


# -- 5. Branch lengths and correlation matrix --------------------------
# Grafen (1989) branch lengths, via ape.

tree_bl <- ape::compute.brlen(tree, method = "Grafen")

phylo_cor_matrix <- ape::vcv(tree_bl, corr = TRUE)


# -- 6. Tip metadata for the figure ------------------------------------

tip_classes <- stats::setNames(sp$Class, sp$Latin.name)
tip_classes <- tip_classes[tree$tip.label]
tip_classes[is.na(tip_classes)] <- "Other"

# Masterdata lists this orchid as Magnoliopsida; orchids are Liliopsida
for (s in c("Cypripedium acaule"))
  if (s %in% names(tip_classes)) tip_classes[s] <- "Liliopsida"


# -- 7. Diagnostics ----------------------------------------------------

missing_from_tree <- setdiff(unique(df$Latin.name[!is.na(df$Latin.name)]),
                             rownames(phylo_cor_matrix))
if (length(missing_from_tree) > 0) {
  cat("\nIn df but NOT in the tree:\n")
  print(missing_from_tree)
} else {
  cat("\nAll species in df have a tip in the tree.\n")
}

cat("phylo_cor_matrix:", nrow(phylo_cor_matrix), "x",
    ncol(phylo_cor_matrix), "\n")
cat("Classes represented:", dplyr::n_distinct(tip_classes), "\n")

saveRDS(list(tree = tree, tree_bl = tree_bl,
             phylo_cor_matrix = phylo_cor_matrix,
             tip_classes = tip_classes,
             phylo_standins = phylo_standins),
        file.path(PATHS$models, "phylogeny.rds"))