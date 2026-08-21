# This script extracts topological N-glycan branching motifs from
# `glydb::glydb_structures()`.

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

expand_sialic_acid_variants <- function(branches) {
  expanded_graphs <- map(as.list(branches), function(graph) {
    monos <- igraph::vertex_attr(graph, "mono")
    sialic_nodes <- which(monos %in% c("Neu5Ac", "Neu5Gc"))
    if (length(sialic_nodes) == 0) {
      return(list(graph))
    }

    assignments <- expand.grid(
      rep(list(c("Neu5Ac", "Neu5Gc")), length(sialic_nodes)),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    map(seq_len(nrow(assignments)), function(id) {
      igraph::set_vertex_attr(
        graph,
        "mono",
        index = sialic_nodes,
        value = as.character(assignments[id, ])
      )
    })
  })
  unique(as_glycan_structure(unlist(expanded_graphs, recursive = FALSE)))
}

topological_glycans <- glydb_structures(
  structure_level = "topological",
  glycan_type = "N"
)
topological_glycans <- topological_glycans[
  !has_floating_parts(topological_glycans) &
    !has_floating_substituents(topological_glycans) &
    !get_alditol(topological_glycans)
]

topological_branches <- extract_branches(topological_glycans)
generic_branches <- as.character(convert_to_generic(topological_branches))
branch_confidence <- attr(topological_branches, "confidence")
best_branch_ids <- vapply(
  split(seq_along(topological_branches), generic_branches),
  function(ids) {
    scores <- branch_confidence[ids]
    scores[is.na(scores)] <- -Inf
    ids[[which.max(scores)]]
  },
  integer(1)
)
best_branch_ids <- sort(unname(best_branch_ids))
topological_branches <- topological_branches[best_branch_ids]
attr(topological_branches, "confidence") <- NULL
topological_branches <- expand_sialic_acid_variants(topological_branches)
topological_branch_index <- split(
  seq_along(topological_branches),
  as.character(convert_to_generic(topological_branches))
)
topological_branch_templates <- map(
  seq_along(topological_branches),
  function(id) {
    branch <- topological_branches[id]
    graph <- as.list(branch)[[1]]
    list(
      nodes = structure_nodes(branch),
      edges = structure_edges(branch),
      root = which(igraph::degree(graph, mode = "in") == 0)
    )
  }
)
names(topological_branch_templates) <- as.character(
  seq_along(topological_branch_templates)
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

default_strucs <- glydb_structures(structure_level = "intact")
default_compositions <- as_glycan_composition(default_strucs)
default_comp_to_struc_metadata <- list(
  structure_keys = as.character(default_strucs),
  floating = has_floating_parts(default_strucs) |
    has_floating_substituents(default_strucs),
  composition = default_compositions,
  generic_keys = as.character(convert_to_generic(default_compositions))
)

usethis::use_data(
  default_comp_to_struc_metadata,
  high_mannose_reference,
  n_glycan_generic_core,
  n_glycan_topological_core,
  topological_branch_index,
  topological_branch_templates,
  topological_branches,
  internal = TRUE,
  overwrite = TRUE
)
