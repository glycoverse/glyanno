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
#'   Defaults to `FALSE`.
#'
#' @returns A tibble with the following columns:
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
  if (is.null(db)) {
    db <- glydb::glydb_compositions(mono_type = "concrete")
  } else {
    db <- .ensure_glycan_composition(db, allow_structure = FALSE)
  }

  # Store confidence before deduplication
  confidences <- attr(db, "confidence")
  db <- unique(db)
  .check_confidence_attr(confidences, return_best)

  # Handle empty composition case
  if (length(comps) == 0) {
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
      confidence = confidences %||% NA_real_
    )
    comps_df <- tibble::tibble(
      composition = glyrepr::convert_to_generic(comps),
      row_id = seq_along(comps)
    )

    res <- comps_df |>
      dplyr::inner_join(db_df, by = c("composition" = "generic")) |>
      dplyr::select(all_of(c(
        "raw" = "composition",
        "enhanced" = "concrete",
        "row_id",
        "confidence"
      )))

    if (return_best) {
      res <- res |>
        dplyr::arrange(.data$row_id, dplyr::desc(.data$confidence)) |>
        dplyr::group_by(.data$row_id) |>
        dplyr::slice(1) |>
        dplyr::ungroup()
    }

    res <- res |>
      dplyr::arrange(.data$row_id) |>
      dplyr::select(all_of(c("raw", "enhanced")))
  } else {
    # Concrete compositions: enhanced = raw
    res <- tibble::tibble(
      raw = comps,
      enhanced = comps
    )
  }

  # Handle zero-row result
  if (nrow(res) == 0) {
    res <- tibble::tibble(
      raw = glyrepr::glycan_composition(),
      enhanced = glyrepr::glycan_composition()
    )
  }
  res
}


#' Check confidence attribute when return_best is TRUE
#'
#' @param confidences The confidence attribute value
#' @param return_best Whether return_best is TRUE
#' @noRd
.check_confidence_attr <- function(confidences, return_best) {
  if (!return_best) {
    return(invisible(NULL))
  }

  if (is.null(confidences)) {
    cli::cli_abort(c(
      "Database must have a {.val confidence} attribute when {.arg return_best} is {.val TRUE}.",
      "i" = "Add confidence scores to the database using {.code attr(db, \"confidence\") <- values}."
    ))
  }

  invisible(NULL)
}
