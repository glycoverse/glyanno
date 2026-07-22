#' Convert glycan composition to glycan structure
#'
#' Given glycan compositions, this function matches them to
#' all possible glycan structures in the `glydb` database.
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
      structure = glyrepr::glycan_structure()
    ))
  }

  mono_type <- glyrepr::get_mono_type(comps)
  db_index <- .prepare_comp_to_struc_index(db, mono_type)
  db <- db_index$db
  .check_return_best_arg(db, return_best)

  if (mono_type == "generic") {
    comp_keys <- glyrepr::convert_to_generic(comps)
  } else {
    if (glyrepr::get_mono_type(db) == "generic") {
      cli::cli_abort(c(
        "Concrete compositions cannot be matched against a generic structure database.",
        "i" = "Use generic compositions (e.g. {.val Hex(1)HexNAc(1)}) or provide a concrete structure database."
      ))
    }
    comp_keys <- comps
  }

  if (return_best) {
    match_ids <- match(
      comp_keys,
      db_index$composition[db_index$best_order]
    )
    return(unname(db[db_index$best_order][match_ids]))
  }

  db_df <- tibble::tibble(
    composition = db_index$composition,
    structure = db,
    confidence = attr(db, "confidence")
  )
  comps_df <- tibble::tibble(
    composition = comp_keys,
    row_id = seq_along(comps)
  )
  res <- comps_df |>
    dplyr::left_join(db_df, by = "composition", relationship = "many-to-many")

  .prepare_result(
    res,
    return_best,
    raw_col = "composition",
    new_col = "structure"
  )
}

.comp_to_struc_cache <- new.env(parent = emptyenv())

.prepare_comp_to_struc_index <- function(db, mono_type) {
  if (!is.null(db)) {
    return(.new_comp_to_struc_index(.prepare_struc_db(db), mono_type))
  }

  cache_key <- paste0("default_", mono_type)
  if (!exists(cache_key, envir = .comp_to_struc_cache, inherits = FALSE)) {
    index <- .new_comp_to_struc_index(
      glydb::glydb_structures(structure_level = "intact"),
      mono_type
    )
    assign(cache_key, index, envir = .comp_to_struc_cache)
  }
  get(cache_key, envir = .comp_to_struc_cache, inherits = FALSE)
}

.new_comp_to_struc_index <- function(db, mono_type) {
  composition <- glyrepr::as_glycan_composition(db)
  if (mono_type == "generic") {
    composition <- glyrepr::convert_to_generic(composition)
  }

  confidence <- attr(db, "confidence")
  best_order <- order(
    replace(confidence, is.na(confidence), -Inf),
    decreasing = TRUE,
    method = "radix"
  )
  list(db = db, composition = composition, best_order = best_order)
}
