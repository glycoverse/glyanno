#' Enhance glycan composition
#'
#' Given a generic or mixed glycan composition (e.g. Hex(5)HexNAc(2)), this
#' function gives all possible compatible concrete glycan compositions (e.g.
#' Man(5)GlcNAc(2)).
#'
#' @inheritSection mz_to_comp How to set `db`
#'
#' @param comps A [glyrepr::glycan_composition()] vector,
#'   or a character vector of glycan composition strings of Byonic or simple style
#'   (e.g. "Hex(5)HexNAc(2)", "H5N4F1S1").
#'   Generic and mixed compositions are matched to compatible concrete
#'   compositions in `db`. Concrete compositions are returned as is.
#' @param db A [glydb::glydb_compositions()] vector,
#'   or a character vector of glycan composition strings of Byonic or simple style
#'   (e.g. "Man(5)GlcNAc(2)", "H5N4F1S1").
#'   All compositions in `db` must be concrete (e.g. Man(5)GlcNAc(2)).
#'   If not provided, `glydb::glydb_compositions(mono_type = "concrete")` will be used.
#' @param return_best Logical. If `TRUE`, only return the highest confidence match
#'   for each input composition. Requires `db` to have a `confidence` attribute.
#'   Use [glydb::glydb_compositions()] for `db` to enable this feature.
#'   Defaults to `FALSE`.
#'
#' @returns If `return_best=TRUE`:
#'   A [glyrepr::glycan_composition()] vector with the same length as `comps`.
#'   Unmatched compositions are returned as `NA`.
#'   If `return_best=FALSE`:
#'   A tibble with the following columns:
#'   - `raw`: The original compositions.
#'   - `enhanced`: The enhanced compositions.
#'   - `confidence`: The database confidence score for each enhanced
#'     composition, or `NA` when no score is available.
#'     Note that one `raw` composition can have different `enhanced` compositions
#'     as multiple rows in the result.
#'
#' @examples
#' enhance_comp("Hex(5)HexNAc(2)")
#'
#' @export
enhance_comp <- function(comps, db = NULL, return_best = FALSE) {
  checkmate::assert_flag(return_best)
  # Input validation and preparation
  comps <- .ensure_glycan_composition(comps, allow_structure = FALSE)
  db_index <- .prepare_enhance_comp_index(db)
  db <- db_index$db
  .check_return_best_arg(db, return_best)

  # Handle empty composition case
  if (length(comps) == 0) {
    if (return_best) {
      return(glyrepr::glycan_composition())
    }
    return(tibble::tibble(
      raw = glyrepr::glycan_composition(),
      enhanced = glyrepr::glycan_composition(),
      confidence = numeric()
    ))
  }

  comps_type <- glyrepr::get_mono_type(comps)
  is_passthrough <- !is.na(comps_type) & comps_type == "concrete"
  is_matchable <- !is_passthrough & !is.na(comps)
  matches <- rep(list(integer()), length(comps))
  if (any(is_matchable)) {
    matchable_comps <- comps[is_matchable]
    matchable_keys <- as.character(matchable_comps)
    unique_ids <- match(matchable_keys, unique(matchable_keys))
    unique_comps <- matchable_comps[!duplicated(matchable_keys)]
    unique_matches <- lapply(
      seq_along(unique_comps),
      function(i) {
        .composition_match_ids(unique_comps[i], db_index$match_index)
      }
    )
    matches[is_matchable] <- unique_matches[unique_ids]
  }

  if (return_best) {
    match_ids <- vapply(
      matches,
      .best_match_id,
      integer(1),
      best_rank = db_index$best_rank
    )
    result <- db[match_ids]
    result[is_passthrough] <- comps[is_passthrough]
    return(unname(result))
  }

  confidence <- attr(db, "confidence") %||% rep(NA_real_, length(db))
  match_lengths <- lengths(matches)
  matched_row_ids <- rep(seq_along(comps), match_lengths)
  candidate_ids <- unlist(matches, use.names = FALSE)
  matched <- tibble::tibble(
    raw = comps[matched_row_ids],
    enhanced = db[candidate_ids],
    confidence = confidence[candidate_ids],
    row_id = matched_row_ids
  )
  passthrough_ids <- which(is_passthrough)
  passthrough <- tibble::tibble(
    raw = comps[passthrough_ids],
    enhanced = comps[passthrough_ids],
    confidence = rep(NA_real_, length(passthrough_ids)),
    row_id = passthrough_ids
  )
  res <- dplyr::bind_rows(matched, passthrough) |>
    dplyr::arrange(.data$row_id)
  .prepare_result(
    res,
    return_best,
    raw_col = "raw",
    new_col = "enhanced"
  )
}

.enhance_comp_cache <- new.env(parent = emptyenv())

.prepare_enhance_comp_index <- function(db) {
  if (!is.null(db)) {
    return(.new_enhance_comp_index(.prepare_comp_db(db)))
  }

  cache_key <- "default"
  if (!exists(cache_key, envir = .enhance_comp_cache, inherits = FALSE)) {
    index <- .new_enhance_comp_index(glydb::glydb_compositions())
    assign(cache_key, index, envir = .enhance_comp_cache)
  }
  get(cache_key, envir = .enhance_comp_cache, inherits = FALSE)
}

.new_enhance_comp_index <- function(db) {
  mono_types <- glyrepr::get_mono_type(db)
  if (any(is.na(mono_types) | mono_types != "concrete")) {
    cli::cli_abort("All compositions in `db` must be concrete.")
  }

  confidence <- attr(db, "confidence")
  best_order <- order(
    replace(confidence, is.na(confidence), -Inf),
    decreasing = TRUE,
    method = "radix"
  )
  best_rank <- match(seq_along(db), best_order)
  list(
    db = db,
    match_index = .new_composition_match_index(db),
    best_rank = best_rank
  )
}
