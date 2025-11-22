#' Convert m/z values to glycan composition
#'
#' Given m/z values, this function matches them to all possible glycan compositions in the `glydb` database.
#'
#' @details
#' # How to set `db`
#'
#' The `db` parameter is very important for all functions in this package.
#' By default, it uses all available glycans in the `glydb` package,
#' which is usually larger than what you need.
#' You can use helper functions in `glydb` to narrow down the database,
#' e.g. [glydb::glydb_compositions()] or [glydb::glydb_structures()].
#'
#' You can use the `species` and `glycan_type` parameters to focus on specific species and glycan type.
#' For example, if you are only interested in N-glycan compositions in human,
#' you can use `glydb::glydb_compositions(species = "Homo sapiens", glycan_type = "N")`.
#' Also, you can decide the level of information in the database by setting `mono_type`
#' of [glydb::glydb_compositions()] and `structure_level` of [glydb::glydb_structures()].
#'
#' You can then pass the result to the `db` parameter of this function.
#' For example,
#'
#' ```
#' my_db <- glydb::glydb_compositions(species = "Homo sapiens", glycan_type = "N")
#' mz_to_comp(mz, db = my_db)
#' ```
#'
#' @param mz A numeric vector of m/z values.
#' @param tol A numeric scalar of the tolerance for the m/z value in Da or a [ppm()] object for dynamic tolerance.
#'   Default is `ppm(10)`.
#' @param db Glycan compositions to match against.
#'   Can be a [glyrepr::glycan_composition()] vector or glycan composition strings
#'   in Byonic style (e.g. Hex(5)HexNAc(2)) or simple style (e.g. H5N4F1S1).
#'   If not provided, `glydb::glydb_compositions(mono_type = "concrete")` will be used.
#' @inheritParams calculate_mz
#'
#' @returns A tibble with the following columns:
#'   - `mz`: The molecule m/z values, same as the input `mz`.
#'   - `composition`: The possible glycan compositions, as [glyrepr::glycan_composition()] vector.
#'   Note that one m/z value can have multiple rows in the result,
#'   corresponding to different possible glycan compositions.
#'
#' @examples
#' mz_to_comp(933.3175, charge = 1, adduct = "Na+")
#'
#' @seealso [ppm()], [glyanno_mass_dict()]
#' @export
mz_to_comp <- function(
  mz,
  tol = ppm(10),
  db = NULL,
  charge = 1,
  adduct = "H+",
  mass_dict = NULL
) {
  checkmate::assert_numeric(mz)
  mz <- mz[!is.na(mz)]
  if (length(mz) == 0) {
    return(tibble::tibble(mz = numeric(0), composition = glyrepr::glycan_composition()))
  }
  checkmate::assert(
    checkmate::check_number(tol),
    checkmate::check_class(tol, "ppm"),
    combine = "or"
  )
  .check_charge_and_adduct(charge, adduct)
  if (!is.null(mass_dict)) {
    .check_custom_mass_dict(mass_dict)
  } else {
    mass_dict <- glyanno_mass_dict(deriv = "none", mass_type = "mono")
  }

  if (is.null(db)) {
    comps <- glydb::glydb_compositions(mono_type = "concrete")
    suppressWarnings(db_mz <- calculate_mz(comps, charge = charge, adduct = adduct, mass_dict = mass_dict, safe = FALSE))
  } else {
    db <- .ensure_glycan_composition(db, allow_structure = FALSE)
    comps <- unique(db)
    suppressWarnings(db_mz <- calculate_mz(comps, charge = charge, adduct = adduct, mass_dict = mass_dict, safe = FALSE))
    na_count <- sum(is.na(db_mz))
    if (na_count > 0) {
      cli::cli_warn(c(
        "Cannot calculate m/z values for {.val {na_count}} glycans in the database.",
        "i" = "They will be dropped before matching."
      ))
    }
  }
  comps <- comps[!is.na(db_mz)]
  db_mz <- db_mz[!is.na(db_mz)]
  db_df <- tibble::tibble(composition = comps, mz = db_mz) |>
    dplyr::mutate(
      tol = ifelse(is.numeric(.env$tol), .env$tol, .env$tol(.data$mz)),
      upper = .data$mz + .data$tol,
      lower = .data$mz - .data$tol
    )

  find_one <- function(mz) {
    db_df |>
      dplyr::filter(.env$mz > .data$lower & .env$mz < .data$upper) |>
      dplyr::pull(.data$composition)
  }

  res_comps <- purrr::map(mz, find_one)
  res_df <- tibble::tibble(mz = mz, composition = res_comps)
  tidyr::unnest(res_df, all_of("composition"))
}
