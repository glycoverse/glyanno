#' Enhance glycan structure
#'
#' Given a glycan structure vector of any resolution level (see
#' [glyrepr::get_structure_level()] for details), this function gives compatible
#' structures with more specific residue identities or linkage information.
#'
#' Input and database vectors may mix generic, concrete, and mixed residues, as
#' well as topological, partial, and intact structures. Each database candidate
#' is matched against each input independently.
#'
#' @param strucs A [glyrepr::glycan_structure()] vector,
#'   or a character vector of glycan structure strings supported by
#'   [glyparse::auto_parse()]. Inputs with unresolved floating parts or
#'   substituents are excluded with a warning.
#' @param db A [glydb::glydb_structures()] vector,
#'   or a character vector of glycan structure strings supported by [glyparse::auto_parse()].
#'   Structures with unresolved floating parts or substituents are excluded
#'   with a warning. The default is [glydb::glydb_structures()] at "intact"
#'   level.
#' @param return_best Logical. If `TRUE`, only return the best matching
#'   structure (highest confidence) for each input structure. `db` must have a
#'   `confidence` attribute. Default is `FALSE`.
#'
#' @returns If `return_best=TRUE`:
#'   An unnamed [glyrepr::glycan_structure()] vector with the same length as `strucs`.
#'   Unmatched structures are returned as `NA`.
#'   If `return_best=FALSE`:
#'   A tibble with the following columns:
#'   - `raw`: The original glycan structures.
#'   - `enhanced`: The enhanced glycan structures.
#'   - `confidence`: The database confidence score for each enhanced
#'     structure, or `NA` when no score is available.
#'     Note that one `raw` glycan structure can have different `enhanced`
#'     structures as multiple rows in the result.
#'
#' @examples
#' # From topological level to intact level
#' db_intact <- c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
#' enhance_struc("Gal(??-?)GalNAc(??-", db = db_intact)
#'
#' # Refine generic residues without changing the structure level
#' db_topo <- "Gal(??-?)GalNAc(??-"
#' enhance_struc("Hex(??-?)HexNAc(??-", db = db_topo)
#'
#' # From partial level to intact level
#' enhance_struc("Gal(b1-?)GalNAc(a1-", db = db_intact)
#'
#' @export
enhance_struc <- function(
  strucs,
  db = NULL,
  return_best = FALSE
) {
  # Input validation and preparation
  strucs <- .ensure_glycan_structure(strucs)
  floating_input <- .replace_floating_structures(strucs)
  strucs <- floating_input$strucs
  rejected_floating <- floating_input$rejected
  checkmate::assert_flag(return_best)

  db <- .prepare_struc_db(db)
  .check_return_best_arg(db, return_best)

  if (length(strucs) == 0) {
    if (return_best) {
      return(glyrepr::glycan_structure())
    }
    return(tibble::tibble(
      raw = glyrepr::glycan_structure(),
      enhanced = glyrepr::glycan_structure(),
      confidence = numeric()
    ))
  }

  if (all(is.na(strucs))) {
    if (return_best) {
      return(strucs)
    }
    keep <- !rejected_floating
    return(tibble::tibble(
      raw = strucs[keep],
      enhanced = strucs[keep],
      confidence = rep(NA_real_, sum(keep))
    ))
  }

  mono_types <- glyrepr::get_mono_type(strucs)
  structure_levels <- glyrepr::get_structure_level(strucs)
  is_passthrough <- !is.na(strucs) &
    mono_types == "concrete" &
    structure_levels == "intact"

  is_missing <- is.na(strucs)
  to_enhance <- strucs[!is_missing & !is_passthrough]
  unique_to_enhance <- unique(to_enhance)

  if (length(db) == 0 || length(unique_to_enhance) == 0) {
    matches <- rep(list(integer()), length(unique_to_enhance))
  } else {
    match_matrix <- glymotif::have_motifs(
      db,
      unique_to_enhance,
      alignments = "whole"
    )
    matches <- lapply(
      seq_along(unique_to_enhance),
      function(i) which(match_matrix[, i])
    )
  }
  unique_ids <- match(
    as.character(strucs),
    as.character(unique_to_enhance)
  )
  confidence <- attr(db, "confidence") %||% rep(NA_real_, length(db))
  best_order <- order(
    replace(confidence, is.na(confidence), -Inf),
    decreasing = TRUE,
    method = "radix"
  )
  best_rank <- match(seq_along(db), best_order)

  if (return_best) {
    result <- purrr::map(seq_along(strucs), function(i) {
      if (is_missing[[i]]) {
        return(glyrepr::glycan_structure(NA))
      }
      if (is_passthrough[[i]]) {
        return(strucs[i])
      }
      best_id <- .best_match_id(matches[[unique_ids[[i]]]], best_rank)
      db[best_id]
    })
    return(unname(.combine_glycan_structures(result)))
  }

  res <- purrr::map_dfr(seq_along(strucs), function(i) {
    if (is_missing[[i]]) {
      return(tibble::tibble(
        raw = strucs[integer()],
        enhanced = db[integer()],
        confidence = numeric(),
        row_id = integer()
      ))
    }
    if (is_passthrough[[i]]) {
      return(tibble::tibble(
        raw = strucs[i],
        enhanced = strucs[i],
        confidence = NA_real_,
        row_id = i
      ))
    }
    db_ids <- matches[[unique_ids[[i]]]]
    tibble::tibble(
      raw = rep(strucs[i], length(db_ids)),
      enhanced = db[db_ids],
      confidence = confidence[db_ids],
      row_id = rep(i, length(db_ids))
    )
  })

  res <- res |>
    dplyr::arrange(.data$row_id) |>
    dplyr::select(all_of(c("raw", "enhanced", "confidence", "row_id")))

  # Handle empty result
  if (nrow(res) == 0) {
    res <- tibble::tibble(
      raw = glyrepr::glycan_structure(),
      enhanced = glyrepr::glycan_structure(),
      confidence = numeric(0),
      row_id = integer(0)
    )
  }

  res |> dplyr::select(-all_of("row_id"))
}

#' Reconstruct topological N-glycan structures de novo
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Reconstructs generic topological N-glycans from their core and branches.
#' Every non-missing input must have only generic residues and no linkage
#' information. Inputs that cannot be reconstructed de novo are matched against
#' a fallback database at "topological" level.
#'
#' Topological de-novo enhancement preserves optional core fucose and
#' bisecting GlcNAc residues. Complex N-glycan branches are matched against the
#' internal branch data. For hybrid N-glycans, the all-Hex arm is assigned as
#' mannose and the HexNAc-bearing arm uses branch matching. High-mannose
#' candidates are constrained to core-aligned subtrees of the Glc3Man9
#' precursor.
#'
#' @param strucs A [glyrepr::glycan_structure()] vector, or a character vector
#'   of glycan structure strings supported by [glyparse::auto_parse()]. Every
#'   non-missing element must be generic and topological. Missing values are
#'   allowed and preserved. Inputs with unresolved floating parts or
#'   substituents are replaced with missing values and produce a warning.
#' @param fallback_db A [glydb::glydb_structures()] vector, or a character
#'   vector of glycan structure strings supported by [glyparse::auto_parse()].
#'   The default is [glydb::glydb_structures()] at "topological" level. Every
#'   non-missing candidate must have concrete residues. Linkage information is
#'   removed before fallback matching. Candidates with unresolved floating
#'   parts or substituents are excluded with a warning. Fallback candidates
#'   require a `confidence` attribute.
#'
#' @returns An unnamed [glyrepr::glycan_structure()] vector with the same length
#'   as `strucs`. Unmatched structures are returned as `NA`.
#'
#' @examples
#' n_generic <- glyrepr::n_glycan_core(linkage = FALSE, mono_type = "generic")
#' enhance_struc_denovo(n_generic)
#'
#' @export
enhance_struc_denovo <- function(strucs, fallback_db = NULL) {
  lifecycle::signal_stage("experimental", "enhance_struc_denovo()")
  strucs <- .ensure_glycan_structure(strucs)
  strucs <- .replace_floating_structures(strucs)$strucs
  .enhance_struc_denovo_topological(strucs, fallback_db)
}

.enhance_struc_denovo_topological <- function(strucs, fallback_db) {
  if (length(strucs) == 0) {
    return(glyrepr::glycan_structure())
  }

  if (all(is.na(strucs))) {
    return(strucs)
  }

  .assert_denovo_strucs(strucs)

  if (!is.null(fallback_db)) {
    fallback_db <- .prepare_denovo_struc_db(fallback_db)
  }

  is_missing <- is.na(strucs)
  unique_strucs <- unique(strucs[!is_missing])
  context <- .denovo_n_glycan_context()
  core_matches <- glymotif::match_motif(
    unique_strucs,
    n_glycan_generic_core,
    alignment = "core"
  )
  unique_candidates <- purrr::map2(
    unique_strucs,
    core_matches,
    function(struc, matches) {
      tryCatch(
        .enhance_n_glycan_topological(struc, matches, context),
        error = function(cnd) glyrepr::glycan_structure()
      )
    }
  )
  unique_candidates <- .materialize_topological_n_glycan_candidates(
    unique_candidates
  )

  unresolved <- lengths(unique_candidates) == 0
  if (any(unresolved)) {
    fallback_db <- if (is.null(fallback_db)) {
      .prepare_denovo_struc_db(NULL)
    } else {
      fallback_db
    }
    .check_return_best_arg(fallback_db, TRUE, arg = "fallback_db")
    fallback_strucs <- unique_strucs[unresolved]
    fallback <- enhance_struc(
      fallback_strucs,
      db = fallback_db,
      return_best = TRUE
    )

    unresolved_ids <- which(unresolved)
    for (i in seq_along(unresolved_ids)) {
      unique_candidates[[unresolved_ids[[i]]]] <- fallback[i]
    }
  }

  unique_keys <- as.character(unique_strucs)
  struc_keys <- as.character(strucs)
  candidates <- purrr::map(seq_along(strucs), function(i) {
    if (is_missing[[i]]) {
      return(glyrepr::glycan_structure(NA))
    }
    unique_candidates[[match(struc_keys[[i]], unique_keys)]]
  })

  best <- purrr::map(candidates, function(x) {
    if (length(x) == 0) {
      return(glyrepr::glycan_structure(NA))
    }
    x[1]
  })
  .combine_glycan_structures(best)
}

.combine_glycan_structures <- function(strucs) {
  iupacs <- unname(purrr::map_chr(strucs, as.character))
  graph_ids <- which(!is.na(iupacs) & !duplicated(iupacs))
  graphs <- purrr::map(
    strucs[graph_ids],
    glyrepr::get_structure_graphs
  )
  names(graphs) <- iupacs[graph_ids]

  glyrepr::new_glycan_structure(iupacs, graphs)
}

.enhance_n_glycan_topological <- function(
  struc,
  core_matches,
  context
) {
  if (length(core_matches) == 0) {
    return(glyrepr::glycan_structure())
  }

  input_graph <- as.list(struc)[[1]]
  core_map <- core_matches[[1]]
  topology <- .classify_n_glycan_topology(
    input_graph,
    core_map,
    context$terminal_core
  )
  core_additions <- .n_glycan_core_additions(
    input_graph,
    core_map,
    context
  )

  if (topology == "high-mannose") {
    return(.enhance_high_mannose_topological(
      struc,
      core_additions
    ))
  }

  branches <- .n_glycan_branch_occurrences(
    input_graph,
    core_map,
    context$terminal_core
  )
  mannose_extensions <- if (topology == "hybrid") {
    .n_glycan_mannose_extensions(
      input_graph,
      core_map,
      context$terminal_core
    )
  } else {
    list()
  }

  if (any(purrr::map_int(branches, \(x) length(x$candidate_ids)) == 0)) {
    return(glyrepr::glycan_structure())
  }

  base <- .topological_n_glycan_base(
    core_additions,
    mannose_extensions,
    context$topological_nodes,
    context$topological_edges
  )
  branch_templates <- topological_branch_templates

  branch_ids <- purrr::map_int(branches, \(x) x$candidate_ids[[1]])
  candidate <- .build_topological_n_glycan(
    branch_ids,
    branches,
    base,
    branch_templates
  )
  .new_topological_n_glycan_candidates(list(candidate))
}

.denovo_n_glycan_context <- function() {
  generic_core_graph <- as.list(n_glycan_generic_core)[[1]]
  out_degree <- igraph::degree(generic_core_graph, mode = "out")
  terminal_core <- which(out_degree == 0)
  central_core <- which(out_degree == 2)
  reducing_core <- .glycan_graph_root(generic_core_graph)

  list(
    terminal_core = terminal_core,
    central_core = central_core,
    reducing_core = reducing_core,
    inner_core = setdiff(
      seq_len(igraph::vcount(generic_core_graph)),
      c(terminal_core, central_core, reducing_core)
    ),
    topological_nodes = glyrepr::structure_nodes(n_glycan_topological_core),
    topological_edges = glyrepr::structure_edges(n_glycan_topological_core)
  )
}

.classify_n_glycan_topology <- function(input_graph, core_map, terminal_core) {
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
  candidate <- igraph::induced_subgraph(reference_graph, matches[[1]]) |>
    glyrepr::glycan_structure() |>
    glyrepr::remove_linkages()
  if (nrow(core_additions) > 0) {
    candidate <- .add_n_glycan_core_additions(candidate, core_additions)
  }
  candidate
}

.high_mannose_reference <- function() {
  high_mannose_reference
}

.n_glycan_mannose_extensions <- function(
  input_graph,
  core_map,
  terminal_core
) {
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

.n_glycan_core_additions <- function(input_graph, core_map, context) {
  extra_children <- function(core_node) {
    children <- igraph::neighbors(
      input_graph,
      core_map[[core_node]],
      mode = "out"
    )
    setdiff(as.integer(children), core_map)
  }

  bisecting <- extra_children(context$central_core)
  core_fuc <- extra_children(context$reducing_core)
  unexpected <- unlist(
    purrr::map(context$inner_core, extra_children),
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
      parent_node = context$central_core,
      input_node = bisecting,
      mono = "GlcNAc",
      sub = ""
    )
  }
  if (length(core_fuc) == 1) {
    .check_core_addition(input_graph, core_fuc, "dHex", "core fucose")
    additions <- tibble::add_row(
      additions,
      parent_node = context$reducing_core,
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
  input_graph,
  core_map,
  terminal_core
) {
  occurrences <- list()

  for (core_node in terminal_core) {
    children <- igraph::neighbors(
      input_graph,
      core_map[[core_node]],
      mode = "out"
    )
    children <- setdiff(as.integer(children), core_map)
    child_monos <- igraph::vertex_attr(input_graph, "mono", index = children)
    branch_roots <- children[child_monos == "HexNAc"]

    for (branch_root in branch_roots) {
      subtree_nodes <- igraph::subcomponent(
        input_graph,
        branch_root,
        mode = "out"
      )
      pattern <- igraph::induced_subgraph(input_graph, subtree_nodes)
      pattern_key <- glyrepr::graph_to_iupac(pattern)
      candidate_ids <- topological_branch_index[[pattern_key]]
      if (is.null(candidate_ids)) {
        candidate_ids <- integer()
      }
      occurrences[[length(occurrences) + 1]] <- list(
        parent_node = core_node,
        candidate_ids = unname(candidate_ids)
      )
    }
  }
  occurrences
}

.topological_n_glycan_base <- function(
  core_additions,
  fixed_extensions,
  core_nodes,
  core_edges
) {
  nodes <- core_nodes
  edges <- core_edges

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
    template <- .topological_subtree_template(extension$structure)
    state <- .append_topological_subtree(
      nodes,
      edges,
      template,
      extension$parent_node
    )
    nodes <- state$nodes
    edges <- state$edges
  }

  list(nodes = nodes, edges = edges)
}

.build_topological_n_glycan <- function(
  branch_ids,
  branches,
  base,
  branch_templates
) {
  nodes <- base$nodes
  edges <- base$edges

  for (i in seq_along(branch_ids)) {
    branch <- branch_templates[[as.character(branch_ids[[i]])]]
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
  list(nodes = nodes, edges = edges)
}

.topological_n_glycan_from_tables <- function(candidates) {
  anomer <- glyrepr::get_anomer(n_glycan_topological_core)
  graphs <- purrr::map(candidates, .topological_n_glycan_graph, anomer = anomer)
  graphs <- purrr::map(graphs, function(graph) {
    graph |>
      glyrepr::validate_glycan_graph() |>
      glyrepr::canonicalize_glycan_graph()
  })
  glyrepr::validate_glycan_graph_vector(graphs)

  iupacs <- purrr::map_chr(graphs, glyrepr::graph_to_iupac)
  keep <- !duplicated(iupacs)
  unique_graphs <- graphs[keep]
  names(unique_graphs) <- iupacs[keep]
  glyrepr::new_glycan_structure(iupacs, unique_graphs)
}

.new_topological_n_glycan_candidates <- function(candidates) {
  structure(candidates, class = "glyanno_topological_n_glycan_candidates")
}

.materialize_topological_n_glycan_candidates <- function(results) {
  candidate_ids <- which(purrr::map_lgl(
    results,
    inherits,
    "glyanno_topological_n_glycan_candidates"
  ))
  if (length(candidate_ids) == 0) {
    return(results)
  }

  counts <- lengths(results[candidate_ids])
  candidates <- unlist(
    purrr::map(results[candidate_ids], unclass),
    recursive = FALSE
  )
  structures <- .topological_n_glycan_from_tables(candidates)
  ends <- cumsum(counts)
  starts <- ends - counts + 1L
  for (i in seq_along(candidate_ids)) {
    results[[candidate_ids[[i]]]] <- structures[starts[[i]]:ends[[i]]]
  }
  results
}

.topological_n_glycan_graph <- function(candidate, anomer) {
  nodes <- candidate$nodes[order(candidate$nodes$node_id), , drop = FALSE]
  edges <- candidate$edges[order(candidate$edges$edge_id), , drop = FALSE]

  graph <- igraph::make_empty_graph(n = nrow(nodes), directed = TRUE)
  if (nrow(edges) > 0) {
    graph <- igraph::add_edges(
      graph,
      c(rbind(edges$from_node, edges$to_node))
    )
  }
  graph <- igraph::set_vertex_attr(
    graph,
    "name",
    value = as.character(seq_len(nrow(nodes)))
  )
  graph <- igraph::set_vertex_attr(graph, "mono", value = nodes$mono)
  graph <- igraph::set_vertex_attr(graph, "sub", value = nodes$sub)
  graph <- igraph::set_edge_attr(graph, "linkage", value = edges$linkage)
  igraph::set_graph_attr(
    graph,
    "anomer",
    anomer
  )
}

.topological_subtree_template <- function(subtree) {
  list(
    nodes = glyrepr::structure_nodes(subtree),
    edges = glyrepr::structure_edges(subtree),
    root = .glycan_graph_root(as.list(subtree)[[1]])
  )
}

.append_topological_subtree <- function(nodes, edges, subtree, parent_node) {
  subtree_nodes <- subtree$nodes
  subtree_edges <- subtree$edges
  subtree_root <- subtree$root
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
