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
#'   Use [glydb::glydb_compositions()] for `db` to enable this feature.
#'   Default is `FALSE`.
#' @inheritParams calculate_mz
#'
#' @returns If `return_best=TRUE`:
#'   An unnamed [glyrepr::glycan_composition()] vector with the same length as `mz`.
#'   Unmatched m/z values are returned as `NA`.
#'   If `return_best=FALSE`:
#'   A tibble with the following columns:
#'   - `mz`: The molecule m/z values, same as the input `mz`.
#'   - `composition`: The possible glycan compositions, as [glyrepr::glycan_composition()] vector.
#'     Note that one m/z value can have multiple rows in the result,
#'     corresponding to different possible glycan compositions.
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
  # Track NA positions so return_best=TRUE can return a same-length vector
  na_input_mask <- is.na(mz)
  mz_to_process <- mz[!na_input_mask]
  if (length(mz_to_process) == 0) {
    if (return_best) {
      return(glyrepr::as_glycan_composition(rep(NA_character_, length(mz))))
    }
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

  db <- .prepare_comp_db(db)
  suppressWarnings(
    db_mz <- calculate_mz(
      db,
      charge = charge,
      adduct = adduct,
      mass_dict = mass_dict,
      safe = FALSE
    )
  )
  na_count <- sum(is.na(db_mz))
  if (na_count > 0 && !.is_glydb_vector(db)) {
    cli::cli_warn(c(
      "Cannot calculate m/z values for {.val {na_count}} glycans in the database.",
      "i" = "They will be dropped before matching."
    ))
  }

  .check_return_best_arg(db, return_best)
  confidences <- attr(db, "confidence")
  # Store the NA mask before filtering
  db_na_mask <- is.na(db_mz)
  db <- db[!db_na_mask]
  db_mz <- db_mz[!db_na_mask]
  # Filter confidences along with db
  if (!is.null(confidences)) {
    confidences <- confidences[!db_na_mask]
  }
  tolerances <- if (is.numeric(tol)) {
    rep(tol, length(db_mz))
  } else {
    tol(db_mz)
  }
  upper <- db_mz + tolerances
  lower <- db_mz - tolerances

  find_one <- function(mz) {
    match_ids <- which(mz > lower & mz < upper)
    if (length(match_ids) == 0) {
      return(integer())
    }
    if (return_best && length(match_ids) > 1) {
      match_confidences <- confidences[match_ids]
      match_confidences[is.na(match_confidences)] <- -Inf
      return(match_ids[which.max(match_confidences)])
    }
    match_ids
  }

  match_ids <- purrr::map(mz_to_process, find_one)
  if (return_best) {
    best_ids <- purrr::map_int(
      match_ids,
      \(ids) if (length(ids) == 0) NA_integer_ else ids[[1]]
    )
    full_ids <- rep(NA_integer_, length(mz))
    full_ids[!na_input_mask] <- best_ids
    unname(db[full_ids])
  } else {
    match_counts <- lengths(match_ids)
    tibble::tibble(
      mz = rep(mz_to_process, match_counts),
      composition = db[unlist(match_ids, use.names = FALSE)]
    )
  }
}
