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
#' @param return_best A logical scalar. If `TRUE`, only the match with the highest confidence
#'   score is returned for each m/z value. The `db` must have a `confidence` attribute.
#'   Default is `FALSE`.
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
  return_best = FALSE,
  charge = 1,
  adduct = "H+",
  mass_dict = NULL
) {
  checkmate::assert_numeric(mz)
  checkmate::assert_flag(return_best)
  mz <- mz[!is.na(mz)]
  if (length(mz) == 0) {
    return(tibble::tibble(
      mz = numeric(0),
      composition = glyrepr::glycan_composition()
    ))
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
    suppressWarnings(
      db_mz <- calculate_mz(
        comps,
        charge = charge,
        adduct = adduct,
        mass_dict = mass_dict,
        safe = FALSE
      )
    )
    # Store confidence BEFORE filtering NA compositions
    confidences <- attr(comps, "confidence")
  } else {
    db <- .ensure_glycan_composition(db, allow_structure = FALSE)
    # Store confidence BEFORE calling unique() (which strips attributes)
    confidences <- attr(db, "confidence")
    comps <- unique(db)
    suppressWarnings(
      db_mz <- calculate_mz(
        comps,
        charge = charge,
        adduct = adduct,
        mass_dict = mass_dict,
        safe = FALSE
      )
    )
    na_count <- sum(is.na(db_mz))
    if (na_count > 0) {
      cli::cli_warn(c(
        "Cannot calculate m/z values for {.val {na_count}} glycans in the database.",
        "i" = "They will be dropped before matching."
      ))
    }
  }
  # Store the NA mask before filtering
  na_mask <- is.na(db_mz)
  # Check confidence attribute after filtering
  .check_confidence_attr(confidences, return_best)
  comps <- comps[!na_mask]
  db_mz <- db_mz[!na_mask]
  # Filter confidences along with comps
  if (!is.null(confidences)) {
    confidences <- confidences[!na_mask]
  }
  db_df <- tibble::tibble(
    composition = comps,
    mz = db_mz,
    confidence = confidences
  ) |>
    dplyr::mutate(
      tol = ifelse(is.numeric(.env$tol), .env$tol, .env$tol(.data$mz)),
      upper = .data$mz + .data$tol,
      lower = .data$mz - .data$tol
    )

  find_one <- function(mz) {
    matches <- db_df |>
      dplyr::filter(.env$mz > .data$lower & .env$mz < .data$upper)
    if (return_best && nrow(matches) > 1) {
      # Arrange by desc(confidence), treating NA as lowest
      matches <- matches |>
        dplyr::mutate(
          conf_sort = ifelse(is.na(.data$confidence), -Inf, .data$confidence)
        ) |>
        dplyr::arrange(dplyr::desc(.data$conf_sort)) |>
        dplyr::select(-"conf_sort")
      matches <- matches[1, ]
    }
    matches |> dplyr::pull(.data$composition)
  }

  res_comps <- purrr::map(mz, find_one)
  res_df <- tibble::tibble(mz = mz, composition = res_comps)
  tidyr::unnest(res_df, all_of("composition"))
}
