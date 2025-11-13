#' Default Mass Dictionary for Glycan Residues
#'
#' A named numeric vector of the mass of each monosaccharide residue.
#' The names are the monosaccharide residues and other necessary ions or molecules,
#' including "Hex", "HexNAc", "dHex", "dHexNAc", "ddHex", "Pen", "HexA", "HexN",
#' "NeuAc", "NeuGc", "Kdn", "Neu", "H+", "H", "H2O", "K+", "Na+", "NH4+", "H-", "Cl-", "HCO3-", and "red_end".
#' "red_end" is the additional mass of the reducing end of the glycan.
#' The values are the masses in Dalton, with 4 decimal places.
#'
#' @param deriv A character scalar of the derivative to use.
#'   Can be "none", "permethyl", or "peracetyl". Default is "none".
#' @param mass_type "mono" for monoisotopic mass, "average" for average mass.
#'   Default is "mono".
#'
#' @returns A named numeric vector.
#' @examples
#' glyanno_mass_dict(deriv = "none", mass_type = "mono")
#' glyanno_mass_dict(deriv = "permethyl", mass_type = "mono")
#' glyanno_mass_dict(deriv = "none", mass_type = "average")
#'
#' @export
glyanno_mass_dict <- function(deriv = "none", mass_type = "mono") {
  checkmate::assert_choice(deriv, c("none", "permethyl", "peracetyl"))
  checkmate::assert_choice(mass_type, c("mono", "average"))

  variable_mass <- switch(paste(deriv, mass_type, sep = "_"),
    "none_mono" = c(
      "Hex" = 162.0528,
      "HexNAc" = 203.0794,
      "dHex" = 146.0579,
      "dHexNAc" = 187.0798,
      "ddHex" = 130.0584,
      "Pen" = 132.0423,
      "HexA" = 176.03209,
      "HexN" = 143.0536,
      "NeuAc" = 291.0954,
      "NeuGc" = 307.0903,
      "Kdn" = 250.0689,
      "Neu" = 249.0802,
      "red_end" = 18.01056
    ),
    "none_average" = c(
      "Hex" = 162.1424,
      "HexNAc" = 203.195,
      "dHex" = 146.143,
      "dHexNAc" = 187.1948,
      "ddHex" = 130.1428,
      "Pen" = 132.1161,
      "HexA" = 176.1259,
      "HexN" = 143.1452,
      "NeuAc" = 291.2579,
      "NeuGc" = 307.2573,
      "Kdn" = 250.2053,
      "Neu" = 249.2188,
      "red_end" = 18.01524
    ),
    "permethyl_mono" = c(
      "Hex" = 204.0998,
      "HexNAc" = 245.1263,
      "dHex" = 174.0892,
      "dHexNAc" = 215.1111,
      "ddHex" = 144.074,
      "Pen" = 160.0736,
      "HexA" = 218.079,
      "HexN" = 217.1268,
      "NeuAc" = 361.1737,
      "NeuGc" = 391.1842,
      "Kdn" = 320.1472,
      "Neu" = 333.1741,
      "red_end" = 46.0419
    ),
    "permethyl_average" = c(
      "Hex" = 204.223,
      "HexNAc" = 245.2756,
      "dHex" = 174.1968,
      "dHexNAc" = 215.2488,
      "ddHex" = 144.1698,
      "Pen" = 160.1699,
      "HexA" = 218.2066,
      "HexN" = 217.2648,
      "NeuAc" = 361.3923,
      "NeuGc" = 391.4186,
      "Kdn" = 320.3397,
      "Neu" = 333.3808,
      "red_end" = 46.069
    ),
    "peracetyl_mono" = c(
      "Hex" = 288.0845,
      "HexNAc" = 287.1005,
      "dHex" = 230.079,
      "dHexNAc" = 271.101,
      "ddHex" = 172.0689,
      "Pen" = 216.0634,
      "HexA" = 260.0532,
      "HexN" = 329.1064,
      "NeuAc" = 417.1271,
      "NeuGc" = 475.1326,
      "Kdn" = 376.1006,
      "Neu" = 501.1436,
      "red_end" = 102.0317
    ),
    "peracetyl_average" = c(
      "Hex" = 288.2542,
      "HexNAc" = 287.2695,
      "dHex" = 230.2176,
      "dHexNAc" = 271.2688,
      "ddHex" = 172.1798,
      "Pen" = 216.1907,
      "HexA" = 260.2005,
      "HexN" = 329.3048,
      "NeuAc" = 417.3698,
      "NeuGc" = 475.4064,
      "Kdn" = 376.3171,
      "Neu" = 501.4408,
      "red_end" = 102.0898
    )
  )

  fixed_mass <- switch(mass_type,
    "mono" = c(
      "H" = 1.00728,
      "H2O" = 18.01056,
      "H+" = 1.00727,
      "K+" = 38.963707,
      "Na+" = 22.989768,
      "NH4+" = 18.033823,
      "H-" = -1.00728,
      "Cl-" = 34.96885271,
      "HCO3-" = 60.98014364
    ),
    "average" = c(
      "H" = 1.00794,
      "H2O" = 18.01524,
      "H+" = 1.00739,
      "K+" = 39.0983,
      "Na+" = 22.998977,
      "NH4+" = 18.0385,
      "H-" = -1.00794,
      "Cl-" = 35.453,
      "HCO3-" = 61.0168
    )
  )

  c(fixed_mass, variable_mass)
}

#' Check the validity of a cumstom mass dictionary
#'
#' Check if a custom mass dictionary is valid.
#' @param mass_dict A named numeric vector of the mass of each monosaccharide residue.
#' @noRd
.check_custom_mass_dict <- function(mass_dict) {
  checkmate::assert_numeric(mass_dict)
  checkmate::assert_named(mass_dict)
  checkmate::assert_set_equal(names(mass_dict), names(glyanno_mass_dict()))
}