#' Enhance glycan structure
#'
#' Given a glycan structure vector of any resolution level (see
#' [glyrepr::get_structure_level()] for details), this function gives all
#' possible glycan structures of higher resolution level.
#'
#' With `method = "db"`, the target resolution level is determined from `db`.
#' When `db` is `NULL`, the default [glydb::glydb_structures()] at "intact"
#' level is used. With `method = "denovo"`, N-glycans are reconstructed from
#' their core and branches without requiring a structure database.
#'
#' Topological de-novo enhancement preserves optional core fucose and
#' bisecting GlcNAc residues. Complex N-glycan branches are matched against the
#' internal branch data. For hybrid N-glycans, the all-Hex arm is assigned as
#' mannose and the HexNAc-bearing arm uses branch matching. High-mannose
#' candidates are constrained to core-aligned subtrees of the Glc3Man9
#' precursor.
#'
#' @param strucs A [glyrepr::glycan_structure()] vector,
#'   or a character vector of glycan structure strings supported by [glyparse::auto_parse()].
#' @param db A [glydb::glydb_structures()] vector,
#'   or a character vector of glycan structure strings supported by [glyparse::auto_parse()].
#'   If not provided, a default structure vector is loaded from [glydb::glydb_structures()]
#'   at "intact" level.
#'   If `db` has a lower or equal resolution level than `strucs`,
#'   the result will be the same as `strucs` (no enhancement).
#' @param return_best Logical. If `TRUE`, only return the best matching
#'   structure (highest confidence) for each input structure. With
#'   `method = "db"`, `db` must have a `confidence` attribute. With
#'   `method = "denovo"`, branch confidence scores are used. Default is
#'   `FALSE`.
#' @param method `r lifecycle::badge("experimental")` Enhancement method.
#'   `"db"` matches complete structures against `db`. `"denovo"` reconstructs
#'   topological N-glycans from their core and branches.
#'
#' @returns If `return_best=TRUE`:
#'   An unnamed [glyrepr::glycan_structure()] vector with the same length as `strucs`.
#'   Unmatched structures are returned as `NA`.
#'   If `return_best=FALSE`:
#'   A tibble with the following columns:
#'   - `raw`: The original glycan structures.
#'   - `enhanced`: The enhanced glycan structures.
#'     Note that one `raw` glycan structure can have different `enhanced` glycan structures
#'     as multiple rows in the result.
#'
#' @examples
#' # From topological level to intact level
#' db_intact <- c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
#' enhance_struc("Gal(??-?)GalNAc(??-", db = db_intact)
#'
#' # From basic level to topological level
#' db_topo <- "Gal(??-?)GalNAc(??-"
#' enhance_struc("Hex(??-?)HexNAc(??-", db = db_topo)
#'
#' # From partial level to intact level
#' enhance_struc("Gal(b1-?)GalNAc(a1-", db = db_intact)
#'
#' # De-novo enhancement of an N-glycan
#' n_basic <- glyrepr::n_glycan_core(linkage = FALSE, mono_type = "generic")
#' enhance_struc(
#'   n_basic,
#'   method = "denovo",
#'   return_best = TRUE
#' )
#'
#' @export
enhance_struc <- function(
  strucs,
  db = NULL,
  return_best = FALSE,
  method = "db"
) {
  lifecycle::signal_stage("experimental", "enhance_struc(method)")

  # Input validation and preparation
  strucs <- .ensure_glycan_structure(strucs)
  checkmate::assert_flag(return_best)
  checkmate::assert_choice(method, c("db", "denovo"))

  if (method == "denovo") {
    return(.enhance_struc_denovo_topological(strucs, return_best))
  }

  db <- .prepare_struc_db(db)
  .check_return_best_arg(db, return_best)

  if (length(strucs) == 0) {
    if (return_best) {
      return(glyrepr::glycan_structure())
    }
    return(tibble::tibble(
      raw = glyrepr::glycan_structure(),
      enhanced = glyrepr::glycan_structure()
    ))
  }

  if (all(is.na(strucs))) {
    if (return_best) {
      return(strucs)
    }
    return(tibble::tibble(
      raw = strucs,
      enhanced = strucs
    ))
  }

  # Define level ranks (higher number means higher resolution)
  level_ranks <- c("basic" = 1, "topological" = 2, "partial" = 3, "intact" = 4)

  # glyrepr::get_structure_level() returns a scalar for the whole vector.
  db_level <- glyrepr::get_structure_level(db)
  target_rank <- level_ranks[[db_level]]
  struc_level <- glyrepr::get_structure_level(strucs)
  struc_rank <- level_ranks[[struc_level]]

  # Create a dataframe to track original IDs
  strucs_df <- tibble::tibble(
    raw = strucs,
    row_id = seq_along(strucs)
  )

  if (struc_rank >= target_rank) {
    res <- strucs_df |>
      dplyr::mutate(enhanced = .data$raw)
  } else {
    to_enhance <- strucs_df$raw
    # Check matches
    # db are targets (glycans), to_enhance are patterns (motifs)
    matches <- glymotif::have_motifs(db, to_enhance, alignments = "whole")

    # For return_best=TRUE: one row per to_enhance, best match or NA
    # For return_best=FALSE: one row per match
    if (return_best) {
      # Build result for each pattern - one best match or NA
      res_list <- purrr::map(seq_along(to_enhance), function(i) {
        col_matches <- which(matches[, i])
        if (length(col_matches) > 0) {
          # Find best match by confidence
          confs <- attr(db, "confidence")[col_matches]
          best_col <- col_matches[which.max(confs)]
          tibble::tibble(
            raw = to_enhance[i],
            enhanced = db[best_col],
            row_id = strucs_df$row_id[i]
          )
        } else {
          # No match - NA
          tibble::tibble(
            raw = to_enhance[i],
            enhanced = glyrepr::glycan_structure(NA),
            row_id = strucs_df$row_id[i]
          )
        }
      })
      res <- dplyr::bind_rows(res_list)
    } else {
      # matches: rows=db, cols=to_enhance
      match_indices <- which(matches, arr.ind = TRUE)

      if (nrow(match_indices) > 0) {
        # match_indices[, "col"] are indices into to_enhance
        # match_indices[, "row"] are indices into db
        matched_row_ids <- strucs_df$row_id[match_indices[, "col"]]
        matched_raw <- strucs_df$raw[match_indices[, "col"]]
        matched_enhanced <- db[match_indices[, "row"]]

        res <- tibble::tibble(
          raw = matched_raw,
          enhanced = matched_enhanced,
          row_id = matched_row_ids
        )
      } else {
        # No matches found
        res <- tibble::tibble(
          raw = to_enhance[integer(0)],
          enhanced = db[integer(0)],
          row_id = integer(0)
        )
      }
    }
  }

  res <- res |>
    dplyr::arrange(.data$row_id) |>
    dplyr::select(all_of(c("raw", "enhanced", "row_id")))

  # Handle empty result
  if (nrow(res) == 0) {
    res <- tibble::tibble(
      raw = glyrepr::glycan_structure(),
      enhanced = glyrepr::glycan_structure(),
      row_id = integer(0)
    )
  }

  if (return_best) {
    # Return vector with same length as input, NA for unmatched
    return(dplyr::pull(res, .data$enhanced))
  }

  res |> dplyr::select(-all_of("row_id"))
}

.enhance_struc_denovo_topological <- function(strucs, return_best) {
  if (length(strucs) == 0) {
    if (return_best) {
      return(glyrepr::glycan_structure())
    }
    return(.empty_enhance_struc_result())
  }

  if (all(is.na(strucs))) {
    if (return_best) {
      return(strucs)
    }
    return(tibble::tibble(raw = strucs, enhanced = strucs))
  }

  if (glyrepr::get_structure_level(strucs) != "basic") {
    if (return_best) {
      return(unname(strucs))
    }
    return(tibble::tibble(raw = strucs, enhanced = strucs))
  }

  candidates <- purrr::map(seq_along(strucs), function(i) {
    if (is.na(strucs)[[i]]) {
      return(glyrepr::glycan_structure())
    }
    .enhance_n_glycan_topological(strucs[i], return_best)
  })

  if (return_best) {
    graphs <- purrr::map(candidates, function(x) {
      if (length(x) == 0) {
        return(NULL)
      }
      as.list(x[1])[[1]]
    })
    return(do.call(glyrepr::glycan_structure, graphs))
  }

  res <- purrr::map2(candidates, seq_along(strucs), function(x, i) {
    if (length(x) == 0) {
      return(NULL)
    }
    tibble::tibble(
      raw = strucs[rep.int(i, length(x))],
      enhanced = x
    )
  }) |>
    dplyr::bind_rows()

  if (nrow(res) == 0) {
    return(.empty_enhance_struc_result())
  }
  res
}

.enhance_n_glycan_topological <- function(struc, return_best) {
  generic_core <- glyrepr::n_glycan_core(
    linkage = FALSE,
    mono_type = "generic"
  )
  core_matches <- glymotif::match_motif(
    struc,
    generic_core,
    alignment = "core"
  )[[1]]
  if (length(core_matches) == 0) {
    cli::cli_abort(
      "De-novo enhancement currently supports only structures with an N-glycan core."
    )
  }

  input_graph <- as.list(struc)[[1]]
  core_graph <- as.list(generic_core)[[1]]
  core_map <- core_matches[[1]]
  topology <- .classify_n_glycan_topology(input_graph, core_graph, core_map)
  core_additions <- .n_glycan_core_additions(
    input_graph,
    core_graph,
    core_map
  )

  if (topology == "high-mannose") {
    return(.enhance_high_mannose_topological(
      struc,
      return_best,
      core_additions
    ))
  }

  branches <- .n_glycan_branch_occurrences(
    struc,
    input_graph,
    core_graph,
    core_map
  )
  mannose_extensions <- if (topology == "hybrid") {
    .n_glycan_mannose_extensions(input_graph, core_graph, core_map)
  } else {
    list()
  }

  if (any(purrr::map_int(branches, \(x) length(x$candidate_ids)) == 0)) {
    return(glyrepr::glycan_structure())
  }

  if (return_best) {
    confidence <- attr(topological_branches, "confidence")
    # Maximizing the total confidence is separable across branch occurrences.
    branch_ids <- purrr::map_int(branches, function(x) {
      scores <- confidence[x$candidate_ids]
      scores[is.na(scores)] <- -Inf
      x$candidate_ids[[which.max(scores)]]
    })
    return(.build_topological_n_glycan(
      branch_ids,
      branches,
      core_additions,
      mannose_extensions
    ))
  }

  if (length(branches) == 0) {
    return(.build_topological_n_glycan(
      integer(),
      branches,
      core_additions,
      mannose_extensions
    ))
  }

  branch_sets <- expand.grid(
    purrr::map(branches, "candidate_ids"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  graphs <- purrr::map(seq_len(nrow(branch_sets)), function(i) {
    candidate <- .build_topological_n_glycan(
      as.integer(branch_sets[i, ]),
      branches,
      core_additions,
      mannose_extensions
    )
    as.list(candidate)[[1]]
  })

  unique(do.call(glyrepr::glycan_structure, graphs))
}

.classify_n_glycan_topology <- function(input_graph, core_graph, core_map) {
  terminal_core <- which(igraph::degree(core_graph, mode = "out") == 0)
  arm_monos <- purrr::map(terminal_core, function(core_node) {
    children <- igraph::neighbors(
      input_graph,
      core_map[[core_node]],
      mode = "out"
    )
    children <- setdiff(as.integer(children), core_map)
    igraph::vertex_attr(input_graph, "mono", index = children)
  })
  root_monos <- unlist(arm_monos, use.names = FALSE)

  if (any(!root_monos %in% c("Hex", "HexNAc"))) {
    cli::cli_abort(
      "N-glycan arm roots must be {.val Hex} or {.val HexNAc} residues."
    )
  }
  if (!any(root_monos == "HexNAc")) {
    return("high-mannose")
  }
  if (!any(root_monos == "Hex")) {
    return("complex")
  }
  mixed_arm <- purrr::map_lgl(arm_monos, function(monos) {
    all(c("Hex", "HexNAc") %in% monos)
  })
  if (any(mixed_arm)) {
    cli::cli_abort(
      "A hybrid N-glycan must have separate all-Hex and HexNAc-bearing arms."
    )
  }
  "hybrid"
}

.enhance_high_mannose_topological <- function(
  struc,
  return_best,
  core_additions
) {
  reference <- .high_mannose_reference()
  motif <- struc
  if (nrow(core_additions) > 0) {
    motif_graph <- as.list(struc)[[1]]
    motif_graph <- igraph::delete_vertices(
      motif_graph,
      core_additions$input_node
    )
    motif <- glyrepr::glycan_structure(motif_graph)
  }
  matches <- glymotif::match_motif(
    reference,
    motif,
    alignment = "core"
  )[[1]]
  if (length(matches) == 0) {
    return(glyrepr::glycan_structure())
  }

  reference_graph <- as.list(reference)[[1]]
  graphs <- purrr::map(matches, function(nodes) {
    igraph::induced_subgraph(reference_graph, nodes)
  })
  candidates <- do.call(glyrepr::glycan_structure, graphs) |>
    glyrepr::reduce_structure_level(to_level = "topological") |>
    unique()
  if (nrow(core_additions) > 0) {
    graphs <- purrr::map(candidates, function(candidate) {
      enhanced <- .add_n_glycan_core_additions(candidate, core_additions)
      as.list(enhanced)[[1]]
    })
    candidates <- do.call(glyrepr::glycan_structure, graphs)
  }

  if (return_best) {
    return(candidates[1])
  }
  candidates
}

.high_mannose_reference <- function() {
  glyparse::auto_parse(paste0(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)",
    "[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
}

.n_glycan_mannose_extensions <- function(input_graph, core_graph, core_map) {
  terminal_core <- which(igraph::degree(core_graph, mode = "out") == 0)
  extensions <- list()

  for (core_node in terminal_core) {
    children <- igraph::neighbors(
      input_graph,
      core_map[[core_node]],
      mode = "out"
    )
    children <- setdiff(as.integer(children), core_map)
    child_monos <- igraph::vertex_attr(input_graph, "mono", index = children)
    hex_roots <- children[child_monos == "Hex"]

    for (root in hex_roots) {
      subtree_nodes <- igraph::subcomponent(input_graph, root, mode = "out")
      subtree_monos <- igraph::vertex_attr(
        input_graph,
        "mono",
        index = subtree_nodes
      )
      if (any(subtree_monos != "Hex")) {
        cli::cli_abort(
          "The all-Hex arm of a hybrid N-glycan can only contain Hex residues."
        )
      }

      subtree <- igraph::induced_subgraph(input_graph, subtree_nodes)
      subtree <- igraph::set_vertex_attr(
        subtree,
        "mono",
        value = rep("Man", igraph::vcount(subtree))
      )
      if (igraph::ecount(subtree) > 0) {
        subtree <- igraph::set_edge_attr(
          subtree,
          "linkage",
          value = rep("??-?", igraph::ecount(subtree))
        )
      }
      subtree <- igraph::set_graph_attr(subtree, "anomer", "??")
      extensions[[length(extensions) + 1]] <- list(
        parent_node = core_node,
        structure = glyrepr::glycan_structure(subtree)
      )
    }
  }
  extensions
}

.n_glycan_core_additions <- function(input_graph, core_graph, core_map) {
  central_core <- which(igraph::degree(core_graph, mode = "out") == 2)
  reducing_core <- .glycan_graph_root(core_graph)
  terminal_core <- which(igraph::degree(core_graph, mode = "out") == 0)

  extra_children <- function(core_node) {
    children <- igraph::neighbors(
      input_graph,
      core_map[[core_node]],
      mode = "out"
    )
    setdiff(as.integer(children), core_map)
  }

  bisecting <- extra_children(central_core)
  core_fuc <- extra_children(reducing_core)
  inner_core <- setdiff(
    seq_along(core_map),
    c(
      central_core,
      reducing_core,
      terminal_core
    )
  )
  unexpected <- unlist(
    purrr::map(inner_core, extra_children),
    use.names = FALSE
  )

  if (length(bisecting) > 1 || length(core_fuc) > 1 || length(unexpected) > 0) {
    cli::cli_abort(
      "The N-glycan core can only have one core fucose and one bisecting GlcNAc."
    )
  }

  additions <- tibble::tibble(
    parent_node = integer(),
    input_node = integer(),
    mono = character(),
    sub = character()
  )
  if (length(bisecting) == 1) {
    .check_core_addition(input_graph, bisecting, "HexNAc", "bisecting GlcNAc")
    additions <- tibble::add_row(
      additions,
      parent_node = central_core,
      input_node = bisecting,
      mono = "GlcNAc",
      sub = ""
    )
  }
  if (length(core_fuc) == 1) {
    .check_core_addition(input_graph, core_fuc, "dHex", "core fucose")
    additions <- tibble::add_row(
      additions,
      parent_node = reducing_core,
      input_node = core_fuc,
      mono = "Fuc",
      sub = ""
    )
  }
  additions
}

.add_n_glycan_core_additions <- function(struc, core_additions) {
  generic_core <- glyrepr::n_glycan_core(
    linkage = FALSE,
    mono_type = "generic"
  )
  core_map <- glymotif::match_motif(
    struc,
    generic_core,
    alignment = "core"
  )[[1]][[1]]
  nodes <- glyrepr::structure_nodes(struc)
  edges <- glyrepr::structure_edges(struc)

  for (i in seq_len(nrow(core_additions))) {
    new_node <- max(nodes$node_id) + 1L
    parent_node <- core_map[[core_additions$parent_node[[i]]]]
    nodes <- dplyr::bind_rows(
      nodes,
      tibble::tibble(
        glycan_id = 1L,
        node_id = new_node,
        mono = core_additions$mono[[i]],
        sub = core_additions$sub[[i]]
      )
    )
    edges <- dplyr::bind_rows(
      edges,
      tibble::tibble(
        glycan_id = 1L,
        edge_id = nrow(edges) + 1L,
        from_node = parent_node,
        to_node = new_node,
        linkage = "??-?"
      )
    )
  }

  glyrepr::structure_from_tibbles(
    nodes,
    edges,
    anomers = glyrepr::get_anomer(struc)
  )
}

.check_core_addition <- function(graph, node, expected_mono, label) {
  mono <- igraph::vertex_attr(graph, "mono", index = node)
  children <- igraph::neighbors(graph, node, mode = "out")
  if (mono != expected_mono || length(children) != 0) {
    cli::cli_abort(
      "The optional {label} must be a terminal {.val {expected_mono}} residue."
    )
  }
}

.n_glycan_branch_occurrences <- function(
  struc,
  input_graph,
  core_graph,
  core_map
) {
  branch_patterns <- glymotif::extract_branch_motif(struc)
  if (length(branch_patterns) == 0) {
    return(list())
  }

  branch_locations <- glymotif::match_motifs(
    struc,
    glymotif::branch_motifs()
  )
  candidate_matches <- glymotif::have_motifs(
    topological_branches,
    branch_patterns,
    alignments = "whole"
  )
  terminal_core <- which(igraph::degree(core_graph, mode = "out") == 0)
  occurrences <- list()

  for (i in seq_along(branch_patterns)) {
    pattern <- branch_patterns[i]
    pattern_name <- as.character(pattern)
    locations <- branch_locations[[pattern_name]][[1]]
    candidate_ids <- which(candidate_matches[, i])
    pattern_root <- .glycan_graph_root(as.list(pattern)[[1]])

    for (location in locations) {
      branch_root <- location[[pattern_root]]
      parent <- as.integer(igraph::neighbors(
        input_graph,
        branch_root,
        mode = "in"
      ))
      parent_core <- match(parent, core_map)
      if (
        length(parent) != 1 ||
          is.na(parent_core) ||
          !parent_core %in% terminal_core
      ) {
        cli::cli_abort(
          "Each N-glycan branch must be attached to a terminal core mannose."
        )
      }
      occurrences[[length(occurrences) + 1]] <- list(
        parent_node = parent_core,
        candidate_ids = candidate_ids
      )
    }
  }
  occurrences
}

.build_topological_n_glycan <- function(
  branch_ids,
  branches,
  core_additions,
  fixed_extensions = list()
) {
  core <- glyrepr::n_glycan_core(linkage = FALSE, mono_type = "concrete")
  nodes <- glyrepr::structure_nodes(core)
  edges <- glyrepr::structure_edges(core)

  for (i in seq_len(nrow(core_additions))) {
    new_node <- max(nodes$node_id) + 1L
    nodes <- dplyr::bind_rows(
      nodes,
      tibble::tibble(
        glycan_id = 1L,
        node_id = new_node,
        mono = core_additions$mono[[i]],
        sub = core_additions$sub[[i]]
      )
    )
    edges <- dplyr::bind_rows(
      edges,
      tibble::tibble(
        glycan_id = 1L,
        edge_id = nrow(edges) + 1L,
        from_node = core_additions$parent_node[[i]],
        to_node = new_node,
        linkage = "??-?"
      )
    )
  }

  for (extension in fixed_extensions) {
    state <- .append_topological_subtree(
      nodes,
      edges,
      extension$structure,
      extension$parent_node
    )
    nodes <- state$nodes
    edges <- state$edges
  }

  for (i in seq_along(branch_ids)) {
    branch <- topological_branches[branch_ids[[i]]]
    state <- .append_topological_subtree(
      nodes,
      edges,
      branch,
      branches[[i]]$parent_node
    )
    nodes <- state$nodes
    edges <- state$edges
  }

  edges$edge_id <- seq_len(nrow(edges))
  glyrepr::structure_from_tibbles(
    nodes,
    edges,
    anomers = glyrepr::get_anomer(core)
  )
}

.append_topological_subtree <- function(nodes, edges, subtree, parent_node) {
  subtree_nodes <- glyrepr::structure_nodes(subtree)
  subtree_edges <- glyrepr::structure_edges(subtree)
  subtree_root <- .glycan_graph_root(as.list(subtree)[[1]])
  node_offset <- max(nodes$node_id)

  subtree_nodes$node_id <- subtree_nodes$node_id + node_offset
  subtree_edges$from_node <- subtree_edges$from_node + node_offset
  subtree_edges$to_node <- subtree_edges$to_node + node_offset
  nodes <- dplyr::bind_rows(nodes, subtree_nodes)
  edges <- dplyr::bind_rows(
    edges,
    subtree_edges,
    tibble::tibble(
      glycan_id = 1L,
      edge_id = nrow(edges) + nrow(subtree_edges) + 1L,
      from_node = parent_node,
      to_node = subtree_root + node_offset,
      linkage = "??-?"
    )
  )
  list(nodes = nodes, edges = edges)
}

.glycan_graph_root <- function(graph) {
  roots <- which(igraph::degree(graph, mode = "in") == 0)
  if (length(roots) != 1) {
    cli::cli_abort(
      "A glycan structure must have exactly one reducing-end root."
    )
  }
  roots
}

.empty_enhance_struc_result <- function() {
  tibble::tibble(
    raw = glyrepr::glycan_structure(),
    enhanced = glyrepr::glycan_structure()
  )
}
