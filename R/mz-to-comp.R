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
#' @param method A character scalar of the method to use. Can be "denovo" or "database". Default is "database".
#' @param db A [glyrepr::glycan_composition()] vector of glycan compositions to match against.
#'   If not provided, `glydb::glydb_compositions(mono_type = "concrete")` will be used.
#'   Only used when `method` is "database".
#' @inheritParams calculate_mz
#'
#' @returns A tibble with the following columns:
#'   - `mass`: The molecule mass.
#'   - `composition`: The possible glycan compositions, as [glyrepr::glycan_composition()] vector.
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
  db = NULL,
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
  checkmate::assert_class(db, "glyrepr_composition", null.ok = TRUE)
  .check_charge_and_adduct(charge, adduct)
  checkmate::assert_choice(method, c("denovo", "database"))
  if (!is.null(mass_dict)) {
    .check_custom_mass_dict(mass_dict)
  } else {
    mass_dict <- glyanno_mass_dict(deriv = "none", mass_type = "mono")
  }

  res <- switch(method,
    "denovo" = .mz_to_comp_denovo(mz, tol, charge, adduct, mass_dict),
    "database" = .mz_to_comp_database(mz, tol, db, charge, adduct, mass_dict)
  )
  res
}

.mz_to_comp_denovo <- function(mz, tol, charge, adduct, mass_dict) {
  stop("Not implemented")
}

.mz_to_comp_database <- function(mz, tol, db, charge, adduct, mass_dict) {
  if (is.null(db)) {
    comps <- glydb::glydb_compositions(mono_type = "concrete")
    suppressWarnings(db_mz <- calculate_mz(comps, charge = charge, adduct = adduct, mass_dict = mass_dict, safe = FALSE))
  } else {
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
    dplyr::filter(db_df, .env$mz > .data$lower & .env$mz < .data$upper)
  }

  purrr::list_rbind(purrr::map(mz, find_one)) |>
    dplyr::select(all_of(c("composition", "mz")))
}