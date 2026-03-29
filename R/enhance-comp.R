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
  db <- .prepare_comp_db(db)
  .check_return_best_arg(db, return_best)

  # Handle empty composition case
  if (length(comps) == 0) {
    if (return_best) {
      return(glyrepr::glycan_composition())
    }
    return(tibble::tibble(
      raw = glyrepr::glycan_composition(),
      enhanced = glyrepr::glycan_composition()
    ))
  }

  # Validate: all compositions in db must be concrete
  if (glyrepr::get_mono_type(db) == "generic") {
    cli::cli_abort("All compositions in `db` must be concrete.")
  }

  # Check if comps are generic or concrete (all same type)
  comps_type <- glyrepr::get_mono_type(comps)

  if (comps_type == "generic") {
    # Generic compositions: match via generic conversion
    db_df <- tibble::tibble(
      generic = glyrepr::convert_to_generic(db),
      concrete = db,
      confidence = attr(db, "confidence") %||% NA_real_
    )
    comps_df <- tibble::tibble(
      composition = glyrepr::convert_to_generic(comps),
      row_id = seq_along(comps)
    )

    res <- comps_df |>
      dplyr::left_join(db_df, by = c("composition" = "generic")) |>
      dplyr::select(all_of(c(
        "raw" = "composition",
        "enhanced" = "concrete",
        "row_id",
        "confidence"
      )))
    res <- .prepare_result(res, return_best, raw_col = "raw", new_col = "enhanced")
  } else {
    # Concrete compositions: enhanced = raw
    if (return_best) {
      return(comps)
    }
    res <- tibble::tibble(
      raw = comps,
      enhanced = comps
    )
  }
  res
}
