test_that("glyanno_mass_dict returns different masses for different parameters", {
  results <- vector("numeric", 6)
  results[[1]] <- glyanno_mass_dict(deriv = "none", mass_type = "mono")[["Hex"]]
  results[[2]] <- glyanno_mass_dict(deriv = "none", mass_type = "average")[["Hex"]]
  results[[3]] <- glyanno_mass_dict(deriv = "permethyl", mass_type = "mono")[["Hex"]]
  results[[4]] <- glyanno_mass_dict(deriv = "permethyl", mass_type = "average")[["Hex"]]
  results[[5]] <- glyanno_mass_dict(deriv = "peracetyl", mass_type = "mono")[["Hex"]]
  results[[6]] <- glyanno_mass_dict(deriv = "peracetyl", mass_type = "average")[["Hex"]]
  expect_equal(length(unique(results)), 6)
})