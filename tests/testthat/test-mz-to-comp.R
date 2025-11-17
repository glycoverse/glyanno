#' Mass dictionary for test
#'
#' This is the integer version of [glyanno_mass_dict()],
#' used for testing purposes.
#'
#' @returns A named numeric vector.
#' @noRd
mass_dict_for_test <- function() {
  c(
    # Fixed masses (mono)
    "H" = 1,
    "H2O" = 18,
    "H+" = 1,
    "K+" = 39,
    "Na+" = 23,
    "NH4+" = 18,
    "H-" = -1,
    "Cl-" = 35,
    "HCO3-" = 61,
    # Variable masses (none_mono)
    "Hex" = 162,
    "HexNAc" = 203,
    "dHex" = 146,
    "dHexNAc" = 187,
    "ddHex" = 130,
    "Pen" = 132,
    "HexA" = 176,
    "HexN" = 143,
    "NeuAc" = 291,
    "NeuGc" = 307,
    "Kdn" = 250,
    "Neu" = 249,
    "S" = 80,
    "P" = 80,
    "red_end" = 18
  )
}

test_that("mz_to_comp works for a simple one m/z case", {
  simple_db <- glyrepr::glycan_composition(
    c(Hex = 1),
    c(Hex = 2),
    c(Hex = 3)
  )
  result <- mz_to_comp(
    mz = 365,
    adduct = "Na+",
    db = simple_db,
    mass_dict = mass_dict_for_test()
  )
  expected <- tibble::tibble(
    mz = 365,
    composition = glyrepr::glycan_composition(c(Hex = 2)),
  )
  expect_equal(result, expected)
})

test_that("mz_to_comp works for a simple two m/z case", {
  simple_db <- glyrepr::glycan_composition(
    c(Glc = 1),
    c(Glc = 1, Gal = 1),
    c(Glc = 2),
    c(Glc = 3)
  )
  result <- mz_to_comp(
    mz = c(203, 365),
    adduct = "Na+",
    db = simple_db,
    mass_dict = mass_dict_for_test()
  )
  expected <- tibble::tibble(
    mz = c(203, 365, 365),
    composition = simple_db[1:3],
  )
  expect_equal(result, expected)
})

test_that("mz_to_comp works for custom numeric tol", {
  simple_db <- glyrepr::glycan_composition(c(Hex = 2))
  mz <- c(365.5, 366.5, 366)
  expect_equal(
    mz_to_comp(mz = mz, tol = 1, db = simple_db, adduct = "Na+", mass_dict = mass_dict_for_test()),
    tibble::tibble(mz = 365.5, composition = glyrepr::glycan_composition(c(Hex = 2)))
  )
})

test_that("mz_to_comp works for custom ppm tol", {
  simple_db <- glyrepr::glycan_composition(c(Hex = 2))
  mz1 <- 365 + 365 / 1e6 * 51  # outsite tol
  mz2 <- 365 + 365 / 1e6 * 49  # within tol
  expect_equal(
    mz_to_comp(mz = c(mz1, mz2), tol = ppm(50), db = simple_db, adduct = "Na+", mass_dict = mass_dict_for_test()),
    tibble::tibble(mz = c(mz2), composition = simple_db)
  )
})

test_that("mz_to_comp handles db with glycans that cannot be calculated m/z values", {
  simple_db <- glyrepr::glycan_composition(c(Glc = 1), c(Mur = 1))
  expect_warning(
    result <- mz_to_comp(mz = 203, db = simple_db, adduct = "Na+", mass_dict = mass_dict_for_test()),
    "Cannot calculate m/z values for 1 glycans in the database."
  )
  expect_equal(result, tibble::tibble(mz = 203, composition = glyrepr::glycan_composition(c(Glc = 1))))
})

test_that("mz_to_comp returns empty tibble with empty mz", {
  expect_equal(
    mz_to_comp(mz = numeric(0), db = glyrepr::glycan_composition(c(Hex = 2)), adduct = "Na+", mass_dict = mass_dict_for_test()),
    tibble::tibble(mz = numeric(0), composition = glyrepr::glycan_composition())
  )
})

test_that("mz_to_comp returns handles NA", {
  expect_equal(
    mz_to_comp(mz = c(365, NA), db = glyrepr::glycan_composition(c(Hex = 2)), adduct = "Na+", mass_dict = mass_dict_for_test()),
    tibble::tibble(mz = 365, composition = glyrepr::glycan_composition(c(Hex = 2)))
  )
})

test_that("mz_to_comp rejects wrong input types", {
  simple_db <- glyrepr::glycan_composition(c(Hex = 2))
  expect_error(mz_to_comp("365", db = simple_db))
  expect_error(mz_to_comp(365, db = glyrepr::n_glycan_core()))
  expect_error(mz_to_comp(365, db = simple_db, adduct = "Na"))
  expect_error(mz_to_comp(365, db = simple_db, tol = "1"))
  expect_error(mz_to_comp(365, db = simple_db, method = "wrong"))
  expect_error(mz_to_comp(365, db = simple_db, mass_dict = c("Hex" = 1)))
})

test_that("mz_to_comp accepts glycan composition strings as db", {
  simple_db <- c("Hex(1)", "Hex(2)", "Hex(3)")
  result <- mz_to_comp(
    mz = 365,
    adduct = "Na+",
    db = simple_db,
    mass_dict = mass_dict_for_test()
  )
  expected <- tibble::tibble(
    mz = 365,
    composition = glyrepr::glycan_composition(c(Hex = 2)),
  )
  expect_equal(result, expected)
})