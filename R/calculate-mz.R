#' Calculate m/z values of glycans
#'
#' This function calculates m/z values from glycans.
#' Different adducts, mass types, and derivatives are supported.
#' Custom mass dictionaries are also supported.
#'
#' @param glycans Glycans to calculate m/z values from.
#'   Supported inputs include [glyrepr::glycan_structure()], [glyrepr::glycan_composition()],
#'   Byonic style composition strings (e.g. Hex(5)HexNAc(2)),
#'   and any structure strings supported by [glyparse::auto_parse()].
#' @param charge Charge to use. Can be 0, 1, 2, 3, -1, -2, -3, etc. 0 means neutral. Default is 1.
#' @param adduct Adduct to use. Can be "H+", "K+", "Na+", "NH4+", "Cl-", "HCO3-".
#'   Default is "H+".
#'   - When `charge` is 0, `adduct` is ignored.
#'   - When `charge` is positive, `adduct` can only be "H+", "K+", "Na+", "NH4+".
#'   - When `charge` is negative, `adduct` can only be "Cl-", "HCO3-".
#' @param mass_dict A named numeric vector of the mass of each monosaccharide residue.
#'   Default is `glyanno_mass_dict(deriv = "none", mass_type = "mono")`.
#'   If a custom mass dictionary is provided,
#'   please make sure the names of the vector are the same as the names in [glyanno_mass_dict()].
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
  mass_dict = NULL
) {
  # ===== Argument processing =====
  comps <- .ensure_glycan_composition(glycans)
  checkmate::assert_int(charge)
  if (charge > 0) {
    if (!checkmate::test_choice(adduct, c("H+", "K+", "Na+", "NH4+"))) {
      cli::cli_abort("When {.arg charge} is positive, {.arg adduct} can only be 'H+', 'K+', 'Na+', or 'NH4+'.")
    }
  } else if (charge < 0) {
    if (!checkmate::test_choice(adduct, c("Cl-", "HCO3-"))) {
      cli::cli_abort("When {.arg charge} is negative, {.arg adduct} can only be 'Cl-' or 'HCO3-'.")
    }
  }
  if (!is.null(mass_dict)) {
    .check_custom_mass_dict(mass_dict)
  } else {
    mass_dict <- glyanno_mass_dict(deriv = "none", mass_type = "mono")
  }

  # ===== m/z calculation =====
  monos <- c("Hex", "HexNAc", "dHex", "dHexNAc", "ddHex", "Pen", "HexA", "HexN", "NeuAc", "NeuGc", "Kdn", "Neu")
  unsupported_monos <- setdiff(glyrepr::available_monosaccharides("generic"), monos)
  counts <- purrr::map(monos, ~ glyrepr::count_mono(comps, .))
  unsupported_counts <- purrr::map(unsupported_monos, ~ glyrepr::count_mono(comps, .))
  has_unsupported_monos <- purrr::map_int(unsupported_counts, sum) > 0
  if (any(has_unsupported_monos)) {
    cli::cli_abort(c(
      "Unsupported monosaccharides found in the glycans.",
      "x" = "Unsupported monosaccharides: {.val {unsupported_monos[has_unsupported_monos]}}",
      "i" = "Supported monosaccharides: {.val {monos}}",
      "i" = "Use {.fn glyrepr::count_mono} to find the invalid glycans."
    ))
  }
  mono_masses <- purrr::map2(monos, counts, ~ mass_dict[.x] * .y)
  mz <- unname(colSums(do.call(rbind, mono_masses)) + mass_dict[adduct] * abs(charge) + mass_dict["red_end"])
  if (charge != 0) {
    mz <- mz / abs(charge)
  }
  mz
}

#' Ensure the glycans is a [glyrepr::glycan_composition()]
#'
#' Convert all supported inputs into a [glyrepr::glycan_composition()].
#' @param glycans Glycans to process.
#' @returns A [glyrepr::glycan_composition()].
#' @noRd
.ensure_glycan_composition <- function(glycans) {
  if (is.character(glycans)) {
    # ===== Character =====
    # Two cases: 1. Byonic composition strings, 2. structure strings
    # First assume composition strings
    tryCatch(
      return(glyrepr::as_glycan_composition(glycans)),
      error = function(e) NULL
    )
    # Then try to parse structure strings
    tryCatch(
      {
        struc <- glyparse::auto_parse(glycans)
        return(glyrepr::as_glycan_composition(struc))
      },
      error = function(e) NULL
    )
    # If both fail, raise an error
    cli::cli_abort("Cannot parse {.arg glycans} as glycan composition or structure strings.")
  } else if (glyrepr::is_glycan_composition(glycans)) {
    # ===== Glycan Composition =====
    glycans
  } else if (glyrepr::is_glycan_structure(glycans)) {
    # ===== Glycan Structure =====
    glyrepr::as_glycan_composition(glycans)
  } else {
    # ===== Other Types =====
    cli::cli_abort(c(
      "{.arg glycans} must be a character vector, a glycan composition, or a glycan structure.",
      "x" = "Got {.cls {class(glycans)}}."
    ))
  }
}