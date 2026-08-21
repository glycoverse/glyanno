#' Convert glycan composition to glycan structure
#'
#' Given glycan compositions, this function matches them to all compatible
#' glycan structures in the `glydb` database. Generic, concrete, and mixed
#' residue identities are matched residue by residue.
#'
#' @inheritSection mz_to_comp How to set `db`
#'
#' @param comps Glycan compositions to match against. Can be either:
#'   - A [glyrepr::glycan_composition()] vector.
#'   - Byonic style composition strings (e.g. Hex(5)HexNAc(2)).
#'   - Simple style composition strings (e.g. H5N4F1S1).
#' @param db Glycan structures to match against.
#'   Can be a [glyrepr::glycan_structure()] vector or any structure strings
#'   supported by [glyparse::auto_parse()].
#'   Structures with unresolved floating parts or substituents are excluded
#'   with a warning.
#'   If not provided, `glydb::glydb_structures(structure_level = "intact")` will be used.
#' @param return_best If `TRUE`, only return the highest confidence match for each
#'   composition. Requires `db` to have a `confidence` attribute.
#'   Use [glydb::glydb_structures()] for `db` to enable this feature.
#'   Default is `FALSE`.
#'
#' @returns If `return_best=TRUE`:
#'   A [glyrepr::glycan_structure()] vector with the same length as `comps`.
#'   Unmatched compositions are returned as `NA`.
#'   If `return_best=FALSE`:
#'   A tibble with the following columns:
#'   - `composition`: The glycan compositions, as [glyrepr::glycan_composition()] vector.
#'   - `structure`: The possible glycan structures, as [glyrepr::glycan_structure()] vector.
#'   - `confidence`: The database confidence score for each structure, or `NA`
#'     when no score is available.
#'     Note that one glycan composition can have multiple rows in the result,
#'     corresponding to different possible glycan structures.
#'
#' @examples
#' comp_to_struc("H5N2")
#'
#' @seealso [glyparse::auto_parse()]
#' @export
comp_to_struc <- function(comps, db = NULL, return_best = FALSE) {
  checkmate::assert_flag(return_best)
  comps <- .ensure_glycan_composition(comps, allow_structure = FALSE)

  if (length(comps) == 0) {
    db <- .prepare_struc_db(db)
    .check_return_best_arg(db, return_best)
    if (return_best) {
      return(glyrepr::glycan_structure())
    }
    return(tibble::tibble(
      composition = glyrepr::glycan_composition(),
      structure = glyrepr::glycan_structure(),
      confidence = numeric()
    ))
  }

  mono_types <- glyrepr::get_mono_type(comps)
  match_mode <- if (all(is.na(mono_types) | mono_types == "concrete")) {
    "exact"
  } else {
    "compatible"
  }
  db_index <- .prepare_comp_to_struc_index(db, match_mode)
  db <- db_index$db
  .check_return_best_arg(db, return_best)

  comp_keys <- as.character(comps)
  unique_ids <- match(comp_keys, unique(comp_keys))
  unique_comps <- comps[!duplicated(comp_keys)]
  unique_matches <- lapply(
    seq_along(unique_comps),
    function(i) {
      .composition_match_ids(unique_comps[i], db_index$match_index)
    }
  )
  matches <- unique_matches[unique_ids]

  if (return_best) {
    match_ids <- vapply(
      matches,
      .best_match_id,
      integer(1),
      best_rank = db_index$best_rank
    )
    return(unname(db[match_ids]))
  }

  confidence <- attr(db, "confidence") %||% rep(NA_real_, length(db))
  match_lengths <- lengths(matches)
  row_ids <- rep(seq_along(comps), match_lengths)
  candidate_ids <- unlist(matches, use.names = FALSE)
  res <- tibble::tibble(
    composition = comps[row_ids],
    structure = db[candidate_ids],
    confidence = confidence[candidate_ids],
    row_id = row_ids
  )

  .prepare_result(
    res,
    return_best,
    raw_col = "composition",
    new_col = "structure"
  )
}

.comp_to_struc_cache <- new.env(parent = emptyenv())

.prepare_comp_to_struc_index <- function(db, match_mode) {
  if (!is.null(db)) {
    return(.new_comp_to_struc_index(
      .prepare_struc_db(db),
      match_mode
    ))
  }

  db <- .prepare_struc_db(NULL)
  composition_key <- "default_composition"
  if (
    !exists(
      composition_key,
      envir = .comp_to_struc_cache,
      inherits = FALSE
    )
  ) {
    composition <- get(
      "intact_composition",
      envir = .default_struc_db_cache,
      inherits = FALSE
    )
    assign(composition_key, composition, envir = .comp_to_struc_cache)
  }

  cache_key <- paste0("default_", match_mode)
  if (!exists(cache_key, envir = .comp_to_struc_cache, inherits = FALSE)) {
    index <- .new_comp_to_struc_index(
      db,
      match_mode,
      get(
        composition_key,
        envir = .comp_to_struc_cache,
        inherits = FALSE
      ),
      get(
        "intact_generic_keys",
        envir = .default_struc_db_cache,
        inherits = FALSE
      )
    )
    assign(cache_key, index, envir = .comp_to_struc_cache)
  }
  get(cache_key, envir = .comp_to_struc_cache, inherits = FALSE)
}

.new_comp_to_struc_index <- function(
  db,
  match_mode,
  composition = glyrepr::as_glycan_composition(db),
  generic_keys = NULL
) {
  confidence <- attr(db, "confidence")
  best_order <- order(
    replace(confidence, is.na(confidence), -Inf),
    decreasing = TRUE,
    method = "radix"
  )
  best_rank <- match(seq_along(db), best_order)
  match_index <- if (match_mode == "exact") {
    .new_exact_composition_match_index(composition)
  } else {
    .new_composition_match_index(composition, generic_keys)
  }
  list(
    db = db,
    match_index = match_index,
    best_rank = best_rank
  )
}
