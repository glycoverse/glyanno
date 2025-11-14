#' Convert m/z values to glycan composition
#'
#' Given m/z values, this function calculates all possible glycan compositions.
#' Two methods are supported: "denovo" and "database".
#' The "denovo" method uses residues masses to calculate all possible glycan compositions.
#' The "database" method matches the m/z value to the glycan composition in `glydb`.
#'
#' @param mz A numeric vector of m/z values.
#' @param tol A numeric scalar of the tolerance for the m/z value in Da or a [ppm()] object for dynamic tolerance.
#'   Default is `ppm(10)`.
#' @param method A character scalar of the method to use. Can be "denovo" or "database".
#'   If "denovo", compositions with generic monosaccharides (e.g. "Hex", "HexNAc") will be returned.
#'   If "database", compositions with concrete monosaccharides (e.g. "Glc", "Gal") will be returned.
#'   Default is "database".
#' @inheritParams calculate_mz
#'
#' @returns A tibble with the following columns:
#'   - `mass`: The molecule mass.
#'   - `composition`: The possible glycan compositions, as [glyrepr::glycan_composition()] objects.
#'   Note that one mass value can have multiple rows in the result,
#'   corresponding to different possible glycan compositions.
#'
#' @examples
#' mz_to_comp(2368.84)
#'
#' @seealso [ppm()], [glyanno_mass_dict()]
#' @export
mz_to_comp <- function(
  mz,
  tol = ppm(10),
  method = "database",
  charge = 1,
  adduct = "H+",
  mass_dict = NULL
) {
  checkmate::assert_numeric(mz)
  checkmate::assert(
    checkmate::check_number(tol),
    checkmate::check_class(tol, "ppm"),
    combine = "or"
  )
  .check_charge_and_adduct(charge, adduct)
  checkmate::assert_choice(method, c("denovo", "database"))
  if (!is.null(mass_dict)) {
    .check_custom_mass_dict(mass_dict)
  } else {
    mass_dict <- glyanno_mass_dict(deriv = "none", mass_type = "mono")
  }

  res <- switch(method,
    "denovo" = .mz_to_comp_denovo(mz, tol, charge, adduct, mass_dict),
    "database" = .mz_to_comp_database(mz, tol, charge, adduct, mass_dict)
  )
  res
}

.mz_to_comp_denovo <- function(mz, tol, charge, adduct, mass_dict) {
  stop("Not implemented")
}

.mz_to_comp_database <- function(mz, tol, charge, adduct, mass_dict) {
  db <- .mz_comp_db(charge, adduct, mass_dict) |>
    dplyr::mutate(
      tol = ifelse(is.numeric(.env$tol), .env$tol, .env$tol(.data$mz)),
      upper = .data$mz + .data$tol,
      lower = .data$mz - .data$tol
    )

  find_one <- function(mz) {
    dplyr::filter(db, .env$mz > .data$lower & .env$mz < .data$upper)
  }

  purrr::list_rbind(purrr::map(mz, find_one)) |>
    dplyr::select(all_of(c("composition", "mz")))
}

#' Database of glycan compositions
#'
#' A refined version of `glydb::fully_determined_glycans`.
#' Used in `mz_to_comp_database()`.
#'
#' @param charge The charge.
#' @param adduct The adduct.
#' @param mass_dict The mass dictionary.
#'
#' @returns A tibble with the following columns:
#'   - `mz`: The m/z value.
#'   - `composition`: The glycan composition.
#'
#' @noRd
.mz_comp_db <- function(charge, adduct, mass_dict) {
  comps <- unique(glydb::fully_determined_glycans$glycan_composition)
  generic_comps <- glyrepr::convert_to_generic(comps)
  suppressWarnings(mz <- calculate_mz(generic_comps, charge = charge, adduct = adduct, mass_dict = mass_dict, safe = FALSE))
  tibble::tibble(composition = comps, mz = mz) |>
    dplyr::filter(!is.na(mz))
}