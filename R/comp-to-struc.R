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
#' @returns A tibble with the following columns:
#'   - `composition`: The glycan compositions, as [glyrepr::glycan_composition()] vector.
#'   - `structure`: The possible glycan structures, as [glyrepr::glycan_structure()] vector.
#'   Note that one glycan composition can have multiple rows in the result,
#'   corresponding to different possible glycan structures.
#'
#' @examples
#' comp_to_struc("H5N2")
#'
#' @seealso [glyparse::auto_parse()]
#' @export
comp_to_struc <- function(comps, db = NULL, return_best = FALSE) {
  checkmate::assert_flag(return_best)
  comps <- .ensure_glycan_composition(comps, allow_structure = FALSE)

  if (is.null(db)) {
    db <- glydb::glydb_structures(structure_level = "intact")
    confidences <- attr(db, "confidence")
  } else {
    confidences <- attr(db, "confidence")
    is_from_glydb <- !is.null(confidences)
    if (!is_from_glydb) {
      db <- .ensure_glycan_structure(db)
      db <- unique(db)
    }
  }

  .check_confidence_attr(confidences, return_best)

  # Handle empty compositions early (before calling get_mono_type)
  if (length(comps) == 0) {
    return(tibble::tibble(
      composition = glyrepr::glycan_composition(),
      structure = glyrepr::glycan_structure()
    ))
  }

  # Get mono type to determine matching strategy
  mono_type <- glyrepr::get_mono_type(comps)

  if (mono_type == "generic") {
    # For generic compositions, convert both to generic for matching
    db_comps_generic <- glyrepr::convert_to_generic(glyrepr::as_glycan_composition(
      db
    ))
    db_df <- tibble::tibble(
      composition = db_comps_generic,
      structure = db,
      confidence = confidences
    )
    comps_df <- tibble::tibble(
      composition = glyrepr::convert_to_generic(comps),
      row_id = seq_along(comps)
    )
    res <- comps_df |>
      dplyr::inner_join(db_df, by = "composition")
    # Filter to best match if requested
    if (return_best) {
      res <- res |>
        dplyr::arrange(.data$row_id, dplyr::desc(.data$confidence)) |>
        dplyr::group_by(.data$row_id) |>
        dplyr::slice(1) |>
        dplyr::ungroup()
    }
    res <- res |>
      dplyr::arrange(.data$row_id) |>
      dplyr::select(all_of(c("composition", "structure")))
  } else {
    # For concrete compositions, match directly to concrete structures only
    # After glyrepr 0.9.0.9000, db must be homogeneous (all generic or all concrete)
    if (glyrepr::get_mono_type(db) == "generic") {
      # Concrete comps cannot match generic structures
      return(tibble::tibble(
        composition = glyrepr::glycan_composition(),
        structure = glyrepr::glycan_structure()
      ))
    }
    db_concrete_df <- tibble::tibble(
      composition = glyrepr::as_glycan_composition(db),
      structure = db,
      confidence = confidences
    )
    comps_df <- tibble::tibble(
      composition = comps,
      row_id = seq_along(comps)
    )
    res <- comps_df |>
      dplyr::inner_join(db_concrete_df, by = "composition")
    # Filter to best match if requested
    if (return_best) {
      res <- res |>
        dplyr::arrange(.data$row_id, dplyr::desc(.data$confidence)) |>
        dplyr::group_by(.data$row_id) |>
        dplyr::slice(1) |>
        dplyr::ungroup()
    }
    res <- res |>
      dplyr::arrange(.data$row_id) |>
      dplyr::select(all_of(c("composition", "structure")))
  }

  # Ensure zero-row result has the expected glycan vector types
  if (nrow(res) == 0) {
    tibble::tibble(
      composition = glyrepr::glycan_composition(),
      structure = glyrepr::glycan_structure()
    )
  } else {
    res
  }
}
