#' Enhance glycan composition
#'
#' Given a generic glycan composition (e.g. Hex(5)HexNAc(2)),
#' this function gives all possible concrete glycan compositions (e.g. Man(5)GlcNAc(2)).
#'
#' @inheritSection mz_to_comp How to set `db`
#'
#' @param comps A [glyrepr::glycan_composition()] vector,
#'   or a character vector of glycan composition strings of Byonic or simple style
#'   (e.g. "Hex(5)HexNAc(2)", "H5N4F1S1").
#'   Generic compositions (e.g. Hex(5)HexNAc(2)) will be matched to all possible concrete compositions in `db`.
#'   Concrete compositions (e.g. Man(5)GlcNAc(2)) will be returned as is.
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

  # Validate: all compositions in db must be concrete
  if (glyrepr::get_mono_type(db) == "generic") {
    cli::cli_abort("All compositions in `db` must be concrete.")
  }

  # Check if comps are generic or concrete (all same type)
  comps_type <- glyrepr::get_mono_type(comps)

  if (comps_type == "generic") {
    comp_keys <- glyrepr::convert_to_generic(comps)
    if (return_best) {
      match_ids <- match(
        comp_keys,
        db_index$generic[db_index$best_order]
      )
      return(unname(db[db_index$best_order][match_ids]))
    }

    db_df <- tibble::tibble(
      generic = db_index$generic,
      concrete = db,
      confidence = attr(db, "confidence") %||% NA_real_
    )
    comps_df <- tibble::tibble(
      composition = comp_keys,
      row_id = seq_along(comps)
    )

    res <- comps_df |>
      dplyr::left_join(
        db_df,
        by = c("composition" = "generic"),
        relationship = "many-to-many"
      ) |>
      dplyr::select(all_of(c(
        "raw" = "composition",
        "enhanced" = "concrete",
        "row_id",
        "confidence"
      )))
    res <- .prepare_result(
      res,
      return_best,
      raw_col = "raw",
      new_col = "enhanced"
    )
  } else {
    # Concrete compositions: enhanced = raw
    if (return_best) {
      return(comps)
    }
    res <- tibble::tibble(
      raw = comps,
      enhanced = comps,
      confidence = NA_real_
    )
  }
  res
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
  confidence <- attr(db, "confidence")
  best_order <- order(
    replace(confidence, is.na(confidence), -Inf),
    decreasing = TRUE,
    method = "radix"
  )
  list(
    db = db,
    generic = glyrepr::convert_to_generic(db),
    best_order = best_order
  )
}
