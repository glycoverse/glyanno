#' Convert glycan composition to glycan structure
#'
#' Given glycan compositions, this function matches them to
#' all possible glycan structures in the `glydb` database.
#'
#' @details
#' # Note about monosaccharide types
#'
#' See [glyrepr::get_mono_type()] for the definition of monosaccharide types.
#' This function is designed to work with glycans with both generic and concrete monosaccharides.
#' It follows the rules:
#' - Generic compositions in `comps` can match both generic and concrete structures in `db`.
#' - Concrete compositions in `comps` can only match concrete structures in `db`.
#'
#' @param comps Glycan compositions to match against. Can be either:
#'   - A [glyrepr::glycan_composition()] vector.
#'   - Byonic style composition strings (e.g. Hex(5)HexNAc(2)).
#'   - Simple style composition strings (e.g. H5N4F1S1).
#' @param db Glycan structures to match against.
#'   Can be a [glyrepr::glycan_structure()] vector or any structure strings
#'   supported by [glyparse::auto_parse()].
#'   If not provided, `glydb::glydb_structures(structure_level = "intact")` will be used.
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
comp_to_struc <- function(comps, db = NULL) {
  comps <- .ensure_glycan_composition(comps, allow_structure = FALSE)
  if (is.null(db)) {
    db <- glydb::glydb_structures(structure_level = "intact")
  } else {
    db <- .ensure_glycan_structure(db)
  }
  db <- unique(db)

  # Below we deal with generic and concrete compositions in `comps` separately.
  comps_df <- tibble::tibble(
    composition = comps,
    mono_type = glyrepr::get_mono_type(comps),
    row_id = seq_along(comps)  # for ordering the result
  )

  # 1. For concrete compositions, we only need to match them to the concrete structures in `db`.
  concrete_db_df <- tibble::tibble(
    composition = glyrepr::as_glycan_composition(db),
    structure = db
  ) |>
    dplyr::filter(glyrepr::get_mono_type(.data$composition) == "concrete")
  concrete_res <- comps_df |>
    dplyr::filter(.data$mono_type == "concrete") |>
    dplyr::inner_join(concrete_db_df, by = "composition")

  # 2. For generic compositions, we first convert compositions in `db` to generic type,
  # then match them to all structures in `db`.
  generic_db_df <- tibble::tibble(
    composition = glyrepr::convert_to_generic(glyrepr::as_glycan_composition(db)),
    structure = db
  )
  generic_res <- comps_df |>
    dplyr::filter(.data$mono_type == "generic") |>
    dplyr::inner_join(generic_db_df, by = "composition")

  # 3. Combine the results from 1 and 2, and arrange by the order of `comps`.
  res <- dplyr::bind_rows(concrete_res, generic_res) |>
    dplyr::arrange(.data$row_id) |>
    dplyr::select(all_of(c("composition", "structure")))

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