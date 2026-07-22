#' Calculate m/z values of glycans
#'
#' This function calculates m/z values from glycans.
#' Different adducts, mass types, and derivatives are supported.
#' Custom mass dictionaries are also supported.
#'
#' @param glycans Glycans to calculate m/z values from. Valid inputs include:
#'   - A [glyrepr::glycan_structure()] vector.
#'   - A [glyrepr::glycan_composition()] vector.
#'   - Byonic style composition strings (e.g. Hex(5)HexNAc(2)).
#'   - Simple style composition strings (e.g. H5N4F1S1).
#'   - Any structure strings supported by [glyparse::auto_parse()].
#' @param charge Charge to use. Can be 0, 1, 2, 3, -1, -2, -3, etc. 0 means neutral. Default is 1.
#' @param adduct Adduct to use. Can be "H+", "K+", "Na+", "NH4+", "H-", "Cl-", "HCO3-".
#'   Default is "H+".
#'   - When `charge` is 0, `adduct` is ignored.
#'   - When `charge` is positive, `adduct` can only be "H+", "K+", "Na+", "NH4+".
#'   - When `charge` is negative, `adduct` can only be "H-", "Cl-", "HCO3-".
#' @param mass_dict A named numeric vector of the mass of each monosaccharide residue.
#'   Default is `glyanno_mass_dict(deriv = "none", mass_type = "mono")`.
#'   If a custom mass dictionary is provided,
#'   please make sure the names of the vector are the same as the names in [glyanno_mass_dict()].
#' @param safe Whether to raise an error when unsupported monosaccharides or substituents are found.
#'   - If `TRUE` (default), an error will be raised.
#'   - If `FALSE`, a warning will be raised and m/z values for invalid glycans will be set to NA.
#'
#' @returns A numeric vector of m/z values.
#' @seealso [glyanno_mass_dict()]
#'
#' @examples
#' library(glyrepr)
#'
#' # Different input types
#' calculate_mz(glycan_composition(c(Gal = 1, GalNAc = 1)), charge = 0)
#' calculate_mz("Gal(1)GalNAc(1)", charge = 0)
#' calculate_mz(as_glycan_structure("Gal(b1-3)GalNAc(a1-"), charge = 0)
#' calculate_mz("Gal(b1-3)GalNAc(a1-", charge = 0)
#'
#' # For common situation in MALDI-TOF MS
#' calculate_mz("Man(5)GlcNAc(2)", charge = 1,adduct = "Na+")
#'
#' # For common situation in ESI MS
#' calculate_mz("Man(5)GlcNAc(2)", charge = 1, adduct = "H+")
#'
#' # Calculate permethylated m/z values
#' calculate_mz("Man(5)GlcNAc(2)", charge = 0, mass_dict = glyanno_mass_dict(deriv = "permethyl"))
#'
#' # Calculate average mass
#' calculate_mz("Man(5)GlcNAc(2)", charge = 0, mass_dict = glyanno_mass_dict(mass_type = "average"))
#'
#' # Vectorization
#' calculate_mz(c("Man(5)GlcNAc(2)", "Gal(1)GalNAc(1)"), charge = 1, adduct = "Na+")
#'
#' @export
calculate_mz <- function(
  glycans,
  charge = 1,
  adduct = "H+",
  mass_dict = NULL,
  safe = TRUE
) {
  # ===== Argument processing =====
  comps <- .ensure_glycan_composition(glycans)
  .check_charge_and_adduct(charge, adduct)
  if (!is.null(mass_dict)) {
    .check_custom_mass_dict(mass_dict)
  } else {
    mass_dict <- glyanno_mass_dict(deriv = "none", mass_type = "mono")
  }

  # ===== m/z calculation =====
  monos <- intersect(
    glyrepr::available_monosaccharides("generic"),
    names(mass_dict)
  )
  subs <- intersect(glyrepr::available_substituents(), names(mass_dict))
  generic_comps <- glyrepr::convert_to_generic(comps)
  counts <- purrr::map(
    c(monos, subs),
    \(component) glyrepr::count_mono(generic_comps, component)
  )

  bad_monos <- setdiff(glyrepr::available_monosaccharides("generic"), monos)
  bad_subs <- setdiff(glyrepr::available_substituents(), subs)
  bad_components <- c(bad_monos, bad_subs)
  bad_component_counts <- purrr::map(
    bad_components,
    \(component) glyrepr::count_mono(generic_comps, component)
  )
  bad_component_totals <- purrr::map_int(
    bad_component_counts,
    \(counts) sum(counts, na.rm = TRUE)
  )
  has_bad_components <- bad_component_totals > 0
  if (safe) {
    if (any(has_bad_components)) {
      cli::cli_abort(c(
        "Unsupported monosaccharides or substituents found in the glycans.",
        "x" = "Unsupported monosaccharides or substituents: {.val {bad_components[has_bad_components]}}",
        "i" = "Supported monosaccharides or substituents: {.val {c(monos, subs)}}",
        "i" = "Use {.fn glyrepr::count_mono} to find the invalid glycans."
      ))
    }
  } else if (any(has_bad_components)) {
    cli::cli_warn(c(
      "Unsupported monosaccharides or substituents found in the glycans.",
      "x" = "Unsupported monosaccharides or substituents: {.val {bad_components[has_bad_components]}}",
      "i" = "Supported monosaccharides or substituents: {.val {c(monos, subs)}}",
      "i" = "Use {.fn glyrepr::count_mono} to find the invalid glycans.",
      "i" = "m/z values for invalid glycans are set to NA."
    ))
  }

  mono_masses <- purrr::map2(c(monos, subs), counts, ~ mass_dict[.x] * .y)
  mz <- unname(
    colSums(do.call(rbind, mono_masses)) +
      mass_dict[adduct] * abs(charge) +
      mass_dict["red_end"]
  )
  if (charge != 0) {
    mz <- mz / abs(charge)
  }
  if (!safe) {
    mz[rowSums(do.call(cbind, bad_component_counts), na.rm = TRUE) > 0] <- NA
  }
  mz
}

.check_charge_and_adduct <- function(charge, adduct) {
  checkmate::assert_int(charge)
  if (charge > 0) {
    if (!checkmate::test_choice(adduct, c("H+", "K+", "Na+", "NH4+"))) {
      cli::cli_abort(
        "When {.arg charge} is positive, {.arg adduct} can only be 'H+', 'K+', 'Na+', or 'NH4+'."
      )
    }
  } else if (charge < 0) {
    if (!checkmate::test_choice(adduct, c("H-", "Cl-", "HCO3-"))) {
      cli::cli_abort(
        "When {.arg charge} is negative, {.arg adduct} can only be 'H-', 'Cl-', or 'HCO3-'."
      )
    }
  }
}
