# This script extracts topological and intact N-glycan branching motifs
# from `glydb::glydb_structures()`.

library(tidyverse)
library(glyrepr)
library(glymotif) # 0.17.3
library(glydb) # 0.6.0

# Extract branches motifs with confidence
extract_branches <- function(glycans) {
  branch_mat <- have_motifs(glycans, branch_motifs())
  branch_data <- branch_mat |>
    as.data.frame() |>
    rownames_to_column("glycan") |>
    as_tibble() |>
    mutate(glycan_confidence = attr(glycans, "confidence"), .after = glycan) |>
    pivot_longer(
      -c(glycan, glycan_confidence),
      names_to = "motif",
      values_to = "has"
    ) |>
    filter(has) |>
    summarise(confidence = max(glycan_confidence), .by = motif)
  branches <- as_glycan_structure(branch_data$motif)
  attr(branches, "confidence") <- branch_data$confidence
  branches
}

intact_glycans <- glydb_structures(
  structure_level = "intact",
  glycan_type = "N"
)
topological_glycans <- glydb_structures(
  structure_level = "topological",
  glycan_type = "N"
)

intact_branches <- extract_branches(intact_glycans)
topological_branches <- extract_branches(topological_glycans)
topological_branch_index <- split(
  seq_along(topological_branches),
  as.character(reduce_structure_level(topological_branches, "basic"))
)

n_glycan_generic_core <- n_glycan_core(
  linkage = FALSE,
  mono_type = "generic"
)
n_glycan_topological_core <- n_glycan_core(
  linkage = FALSE,
  mono_type = "concrete"
)
high_mannose_reference <- glyparse::auto_parse(paste0(
  "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)",
  "[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]",
  "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
))

usethis::use_data(
  high_mannose_reference,
  intact_branches,
  n_glycan_generic_core,
  n_glycan_topological_core,
  topological_branch_index,
  topological_branches,
  internal = TRUE,
  overwrite = TRUE
)
