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

test_that("mz_to_comp preserves duplicated m/z values", {
  db <- glyrepr::glycan_composition(
    c(Glc = 2),
    c(Gal = 2),
    c(Glc = 1)
  )
  attr(db, "confidence") <- c(1, 2, 3)
  mz <- c(365, 203, 365)

  best <- mz_to_comp(
    mz,
    db = db,
    adduct = "Na+",
    mass_dict = mass_dict_for_test(),
    return_best = TRUE
  )
  expanded <- mz_to_comp(
    mz,
    db = db,
    adduct = "Na+",
    mass_dict = mass_dict_for_test(),
    return_best = FALSE
  )

  expect_equal(as.character(best), c("Gal(2)", "Glc(1)", "Gal(2)"))
  expect_equal(expanded$mz, c(365, 365, 203, 365, 365))
  expect_equal(
    as.character(expanded$composition),
    c("Glc(2)", "Gal(2)", "Glc(1)", "Glc(2)", "Gal(2)")
  )
})

test_that("mz_to_comp works for custom numeric tol", {
  simple_db <- glyrepr::glycan_composition(c(Hex = 2))
  mz <- c(365.5, 366.5, 366)
  expect_equal(
    mz_to_comp(
      mz = mz,
      tol = 1,
      db = simple_db,
      adduct = "Na+",
      mass_dict = mass_dict_for_test()
    ),
    tibble::tibble(
      mz = 365.5,
      composition = glyrepr::glycan_composition(c(Hex = 2))
    )
  )
})

test_that("mz_to_comp works for custom ppm tol", {
  simple_db <- glyrepr::glycan_composition(c(Hex = 2))
  mz1 <- 365 + 365 / 1e6 * 51 # outsite tol
  mz2 <- 365 + 365 / 1e6 * 49 # within tol
  expect_equal(
    mz_to_comp(
      mz = c(mz1, mz2),
      tol = ppm(50),
      db = simple_db,
      adduct = "Na+",
      mass_dict = mass_dict_for_test()
    ),
    tibble::tibble(mz = c(mz2), composition = simple_db)
  )
})

test_that("mz_to_comp handles db with glycans that cannot be calculated m/z values", {
  simple_db <- glyrepr::glycan_composition(c(Glc = 1), c(Mur = 1))
  expect_warning(
    result <- mz_to_comp(
      mz = 203,
      db = simple_db,
      adduct = "Na+",
      mass_dict = mass_dict_for_test()
    ),
    "Cannot calculate m/z values for 1 glycans in the database."
  )
  expect_equal(
    result,
    tibble::tibble(
      mz = 203,
      composition = glyrepr::glycan_composition(c(Glc = 1))
    )
  )
})

test_that("mz_to_comp returns empty tibble with empty mz", {
  expect_equal(
    mz_to_comp(
      mz = numeric(0),
      db = glyrepr::glycan_composition(c(Hex = 2)),
      adduct = "Na+",
      mass_dict = mass_dict_for_test()
    ),
    tibble::tibble(mz = numeric(0), composition = glyrepr::glycan_composition())
  )
})

test_that("mz_to_comp returns handles NA", {
  expect_equal(
    mz_to_comp(
      mz = c(365, NA),
      db = glyrepr::glycan_composition(c(Hex = 2)),
      adduct = "Na+",
      mass_dict = mass_dict_for_test()
    ),
    tibble::tibble(
      mz = 365,
      composition = glyrepr::glycan_composition(c(Hex = 2))
    )
  )
})

test_that("mz_to_comp with return_best=TRUE preserves NA positions in input", {
  db <- glyrepr::glycan_composition(c(Hex = 2))
  attr(db, "confidence") <- c(1.0)

  result <- mz_to_comp(
    mz = c(365, NA, 999),
    adduct = "Na+",
    db = db,
    mass_dict = mass_dict_for_test(),
    return_best = TRUE
  )

  expect_false(is.data.frame(result))
  expect_equal(length(result), 3)
  expect_equal(as.character(result[1]), "Hex(2)")
  expect_true(is.na(result[2]))
  expect_true(is.na(result[3]))
})

test_that("mz_to_comp with return_best=TRUE returns all-NA vector for all-NA input", {
  db <- glyrepr::glycan_composition(c(Hex = 2))
  attr(db, "confidence") <- c(1.0)

  result <- mz_to_comp(
    mz = c(NA_real_, NA_real_),
    adduct = "Na+",
    db = db,
    mass_dict = mass_dict_for_test(),
    return_best = TRUE
  )

  expect_false(is.data.frame(result))
  expect_equal(length(result), 2)
  expect_true(all(is.na(result)))
})

test_that("mz_to_comp rejects wrong input types", {
  simple_db <- glyrepr::glycan_composition(c(Hex = 2))
  expect_error(mz_to_comp("365", db = simple_db))
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

test_that("mz_to_comp with return_best=TRUE keeps highest confidence match", {
  db <- glyrepr::glycan_composition(
    c(Glc = 2),
    c(Gal = 2)
  )
  attr(db, "confidence") <- c(1.0, 2.0)

  result <- mz_to_comp(
    mz = 365,
    adduct = "Na+",
    db = db,
    mass_dict = mass_dict_for_test(),
    return_best = TRUE
  )

  expect_false(is.data.frame(result))
  expect_equal(length(result), 1)
  expect_equal(as.character(result[1]), "Gal(2)")
})

test_that("mz_to_comp with return_best=TRUE returns vector with NA for no match", {
  db <- glyrepr::glycan_composition(c(Glc = 1))
  attr(db, "confidence") <- c(1.0)
  suppressWarnings(
    db_mz <- calculate_mz(
      db,
      charge = 1,
      adduct = "H+",
      mass_dict = mass_dict_for_test()
    )
  )

  # First mz matches Glc(1), second doesn't match anything
  mz_values <- c(db_mz[1], 999.0)
  result <- mz_to_comp(
    mz_values,
    charge = 1,
    adduct = "H+",
    db = db,
    mass_dict = mass_dict_for_test(),
    return_best = TRUE
  )

  expect_false(is.data.frame(result))
  expect_equal(length(result), 2)
  expect_equal(as.character(result[1]), "Glc(1)")
  expect_true(is.na(result[2]))
})

test_that("mz_to_comp with return_best=TRUE errors without confidence attr", {
  db <- glyrepr::glycan_composition(c(Hex = 2))

  expect_error(
    mz_to_comp(
      mz = 365,
      db = db,
      adduct = "Na+",
      mass_dict = mass_dict_for_test(),
      return_best = TRUE
    ),
    "must have a .*confidence.* attribute"
  )
})

test_that("mz_to_comp treats NA confidence as lowest", {
  db <- glyrepr::glycan_composition(
    c(Glc = 2),
    c(Gal = 2)
  )
  attr(db, "confidence") <- c(NA_real_, 2.0)

  result <- mz_to_comp(
    mz = 365,
    adduct = "Na+",
    db = db,
    mass_dict = mass_dict_for_test(),
    return_best = TRUE
  )

  expect_false(is.data.frame(result))
  expect_equal(length(result), 1)
  expect_equal(as.character(result[1]), "Gal(2)")
})

test_that("mz_to_comp works with db=NULL (regression: confidences undefined)", {
  # This test ensures confidences variable is defined when db=NULL
  # Should not error even without return_best
  expect_no_error(
    result <- mz_to_comp(
      mz = 933.3175,
      adduct = "Na+",
      db = NULL,
      return_best = FALSE
    )
  )
  # Result should be a tibble with mz and composition columns
  expect_named(result, c("mz", "composition"))
  # Should have at least one match from default glydb
  expect_gt(nrow(result), 0)
})
