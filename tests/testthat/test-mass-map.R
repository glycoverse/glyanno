test_that("glyanno_mass_dict returns different masses for different parameters", {
  results <- vector("numeric", 6)
  results[[1]] <- glyanno_mass_dict(deriv = "none", mass_type = "mono")[["Hex"]]
  results[[2]] <- glyanno_mass_dict(deriv = "none", mass_type = "average")[[
    "Hex"
  ]]
  results[[3]] <- glyanno_mass_dict(deriv = "permethyl", mass_type = "mono")[[
    "Hex"
  ]]
  results[[4]] <- glyanno_mass_dict(
    deriv = "permethyl",
    mass_type = "average"
  )[["Hex"]]
  results[[5]] <- glyanno_mass_dict(deriv = "peracetyl", mass_type = "mono")[[
    "Hex"
  ]]
  results[[6]] <- glyanno_mass_dict(
    deriv = "peracetyl",
    mass_type = "average"
  )[["Hex"]]
  expect_equal(length(unique(results)), 6)
})

test_that("underivatized HexN uses the hexosamine residue mass", {
  # Unimod accession 454: C6H11NO4, monoisotopic 161.068808,
  # average 161.1558 Da. The dictionary uses four decimal places.
  mono <- glyanno_mass_dict(deriv = "none", mass_type = "mono")
  average <- glyanno_mass_dict(deriv = "none", mass_type = "average")

  expect_equal(unname(mono["HexN"]), 161.0688)
  expect_equal(unname(average["HexN"]), 161.1558)
  expect_equal(calculate_mz("HexN(1)", charge = 0), 179.07936)
  expect_equal(
    calculate_mz("HexN(1)", charge = 0, mass_dict = average),
    179.17104
  )
  expect_equal(
    calculate_mz("HexN(2)Hex(1)", charge = 2),
    (2 * 161.0688 + 162.0528 + 18.01056 + 2 * 1.00727) / 2
  )
})


test_that("underivatized residue masses agree with elemental compositions", {
  # Residue formulas exclude one water molecule from each free sugar.
  # HUPO-PSI ProForma: HexN C6H11NO4 and Neu C9H15NO7.
  # PubChem: dHexNAc free sugar C8H15NO5; ddHex free sugar C6H12O4.
  atoms <- c(
    C = 12,
    H = 1.00782503223,
    N = 14.00307400443,
    O = 15.99491461957
  )
  formulas <- list(
    HexN = c(C = 6, H = 11, N = 1, O = 4),
    dHexNAc = c(C = 8, H = 13, N = 1, O = 4),
    ddHex = c(C = 6, H = 10, O = 3),
    Neu = c(C = 9, H = 15, N = 1, O = 7)
  )
  exact <- vapply(formulas, function(formula) {
    sum(formula * atoms[names(formula)])
  }, numeric(1))
  mass_dict <- glyanno_mass_dict()
  expect_equal(mass_dict[names(formulas)], round(exact, 4))
  for (residue in names(formulas)) {
    expect_equal(
      calculate_mz(paste0(residue, "(1)"), charge = 0),
      unname(round(exact[residue], 4) + mass_dict["red_end"])
    )
  }
})
