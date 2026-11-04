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
