#' Default Mass Dictionary for Glycan Residues
#'
#' A named numeric vector of the mass of each monosaccharide residue.
#' The names are the monosaccharide residues,
#' including "Hex", "HexNAc", "dHex", "Pent", "NeuAc", "NeuGc", "HexA".
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

  switch(paste(deriv, mass_type, sep = "_"),
    "none_mono" = c(
      "Hex" = 162.0528,
      "HexNAc" = 203.0794,
      "dHex" = 146.0579,
      "Pent" = 132.0423,
      "NeuAc" = 291.0954,
      "NeuGc" = 307.0903,
      "HexA" = 176.03209
    ),
    "none_average" = c(
      "Hex" = 162.1424,
      "HexNAc" = 203.195,
      "dHex" = 146.143,
      "Pent" = 132.1161,
      "NeuAc" = 291.2579,
      "NeuGc" = 307.2573,
      "HexA" = 176.1259
    ),
    "permethyl_mono" = c(
      "Hex" = 204.0998,
      "HexNAc" = 245.1263,
      "dHex" = 174.0892,
      "Pent" = 160.0736,
      "NeuAc" = 361.1737,
      "NeuGc" = 391.1842,
      "HexA" = 218.079
    ),
    "permethyl_average" = c(
      "Hex" = 204.223,
      "HexNAc" = 245.2756,
      "dHex" = 174.1968,
      "Pent" = 160.1699,
      "NeuAc" = 361.3923,
      "NeuGc" = 391.4186,
      "HexA" = 218.2066
    ),
    "peracetyl_mono" = c(
      "Hex" = 288.0845,
      "HexNAc" = 287.1005,
      "dHex" = 230.079,
      "Pent" = 216.0634,
      "NeuAc" = 417.1271,
      "NeuGc" = 475.1326,
      "HexA" = 260.0532
    ),
    "peracetyl_average" = c(
      "Hex" = 288.2542,
      "HexNAc" = 287.2695,
      "dHex" = 230.2176,
      "Pent" = 216.1907,
      "NeuAc" = 417.3698,
      "NeuGc" = 475.4064,
      "HexA" = 260.2005
    )
  )
}