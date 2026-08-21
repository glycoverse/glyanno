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

  db_index <- .prepare_comp_to_struc_index(db)
  db <- db_index$db
  .check_return_best_arg(db, return_best)

  if (return_best) {
    match_ids <- vapply(
      seq_along(comps),
      function(i) {
        candidate_ids <- .composition_match_ids(
          comps[i],
          db_index$match_index
        )
        .best_match_id(candidate_ids, db_index$best_rank)
      },
      integer(1)
    )
    return(unname(db[match_ids]))
  }

  confidence <- attr(db, "confidence") %||% rep(NA_real_, length(db))
  res <- purrr::map_dfr(seq_along(comps), function(i) {
    candidate_ids <- .composition_match_ids(comps[i], db_index$match_index)
    tibble::tibble(
      composition = rep(comps[i], length(candidate_ids)),
      structure = db[candidate_ids],
      confidence = confidence[candidate_ids],
      row_id = rep(i, length(candidate_ids))
    )
  })

  .prepare_result(
    res,
    return_best,
    raw_col = "composition",
    new_col = "structure"
  )
}

.comp_to_struc_cache <- new.env(parent = emptyenv())

.prepare_comp_to_struc_index <- function(db) {
  if (!is.null(db)) {
    return(.new_comp_to_struc_index(.prepare_struc_db(db)))
  }

  cache_key <- "default"
  if (!exists(cache_key, envir = .comp_to_struc_cache, inherits = FALSE)) {
    index <- .new_comp_to_struc_index(
      .prepare_struc_db(NULL)
    )
    assign(cache_key, index, envir = .comp_to_struc_cache)
  }
  get(cache_key, envir = .comp_to_struc_cache, inherits = FALSE)
}

.new_comp_to_struc_index <- function(db) {
  composition <- glyrepr::as_glycan_composition(db)
  confidence <- attr(db, "confidence")
  best_order <- order(
    replace(confidence, is.na(confidence), -Inf),
    decreasing = TRUE,
    method = "radix"
  )
  best_rank <- match(seq_along(db), best_order)
  list(
    db = db,
    match_index = .new_composition_match_index(composition),
    best_rank = best_rank
  )
}
