test_that("comp_to_struc works for generic db and generic comps", {
  db <- glyrepr::as_glycan_structure(c(
    "Hex(??-?)[HexNAc(??-?)]HexNAc(??-",
    "Hex(??-?)HexNAc(??-",
    "HexNAc(??-"
  ))
  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  result <- comp_to_struc(comps, db) |>
    dplyr::mutate(structure = as.character(structure))
  expected <- tibble::tibble(
    composition = comps,
    structure = "Hex(??-?)HexNAc(??-"
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc works for concrete db and generic comps", {
  db <- glyrepr::as_glycan_structure(c(
    "GalNAc(a1-",
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  result <- comp_to_struc(comps, db) |>
    dplyr::mutate(structure = as.character(structure))
  expected <- tibble::tibble(
    composition = comps,
    structure = c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc errors for generic db and concrete comps", {
  db <- glyrepr::as_glycan_structure("HexNAc(??-")
  comps <- glyrepr::as_glycan_composition("GalNAc(1)")
  expect_error(
    comp_to_struc(comps, db),
    "Concrete compositions cannot be matched against a generic structure database"
  )
})

test_that("comp_to_struc works for concrete db and concrete comps", {
  db <- glyrepr::as_glycan_structure(c(
    "GalNAc(a1-",
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  comps <- glyrepr::as_glycan_composition("Gal(1)GalNAc(1)")
  result <- comp_to_struc(comps, db) |>
    dplyr::mutate(structure = as.character(structure))
  expected <- tibble::tibble(
    composition = comps,
    structure = c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc works with multiple compositions", {
  db <- glyrepr::as_glycan_structure(c(
    "Hex(??-?)[HexNAc(??-?)]HexNAc(??-",
    "Hex(??-?)HexNAc(??-",
    "HexNAc(??-"
  ))
  comps <- glyrepr::as_glycan_composition(c("Hex(1)HexNAc(1)", "HexNAc(1)"))
  result <- comp_to_struc(comps, db) |>
    dplyr::mutate(structure = as.character(structure))
  expected <- tibble::tibble(
    composition = comps,
    structure = c("Hex(??-?)HexNAc(??-", "HexNAc(??-")
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc works with different composition types", {
  db <- glyrepr::as_glycan_structure(c(
    "Hex(??-?)[HexNAc(??-?)]HexNAc(??-",
    "Hex(??-?)HexNAc(??-",
    "HexNAc(??-"
  ))
  comps <- c("Hex(1)HexNAc(1)", "N1")
  result <- comp_to_struc(comps, db) |>
    dplyr::mutate(structure = as.character(structure))
  expected <- tibble::tibble(
    composition = glyrepr::as_glycan_composition(comps),
    structure = c("Hex(??-?)HexNAc(??-", "HexNAc(??-")
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc handles duplicate compositions", {
  db <- glyrepr::as_glycan_structure(c(
    "Hex(??-?)HexNAc(??-",
    "HexNAc(??-"
  ))
  comps <- glyrepr::as_glycan_composition(c(
    "Hex(1)HexNAc(1)",
    "HexNAc(1)",
    "Hex(1)HexNAc(1)"
  ))
  result <- comp_to_struc(comps, db) |>
    dplyr::mutate(structure = as.character(structure))
  expected <- tibble::tibble(
    composition = comps,
    structure = c("Hex(??-?)HexNAc(??-", "HexNAc(??-", "Hex(??-?)HexNAc(??-")
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc accepts empty compositions", {
  db <- glyrepr::as_glycan_structure(c("Hex(??-?)HexNAc(??-", "HexNAc(??-"))
  comps <- glyrepr::glycan_composition()
  result <- comp_to_struc(comps, db)
  expected <- tibble::tibble(
    composition = glyrepr::glycan_composition(),
    structure = glyrepr::glycan_structure()
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc accepts empty db", {
  db <- glyrepr::glycan_structure()
  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  result <- comp_to_struc(comps, db)
  expected <- tibble::tibble(
    composition = glyrepr::glycan_composition(),
    structure = glyrepr::glycan_structure()
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc with return_best=TRUE keeps highest confidence match", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  attr(db, "confidence") <- c(1.0, 2.0)

  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  result <- comp_to_struc(comps, db, return_best = TRUE)

  expect_equal(as.character(result), "Gal(b1-4)GalNAc(a1-")
})

test_that("comp_to_struc with return_best=TRUE errors without confidence attr", {
  db <- glyrepr::as_glycan_structure("HexNAc(??-")
  comps <- glyrepr::as_glycan_composition("HexNAc(1)")

  expect_error(
    comp_to_struc(comps, db, return_best = TRUE),
    "must have a .*confidence.* attribute"
  )
})

test_that("comp_to_struc with return_best=FALSE keeps all matches (backward compat)", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  attr(db, "confidence") <- c(1.0, 2.0)

  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  result <- comp_to_struc(comps, db, return_best = FALSE) |>
    dplyr::mutate(structure = as.character(structure))

  expect_equal(nrow(result), 2)
})

test_that("comp_to_struc tie-breaking keeps first when confidence equal", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  attr(db, "confidence") <- c(2.0, 2.0)

  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  result <- comp_to_struc(comps, db, return_best = TRUE)

  expect_equal(as.character(result), "Gal(b1-3)GalNAc(a1-")
})

test_that("comp_to_struc treats NA confidence as lowest", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  attr(db, "confidence") <- c(NA_real_, 2.0)

  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  result <- comp_to_struc(comps, db, return_best = TRUE)

  expect_equal(as.character(result), "Gal(b1-4)GalNAc(a1-")
})

test_that("comp_to_struc works with db=NULL (regression: confidences undefined)", {
  # This test ensures confidences variable is defined when db=NULL
  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  # Should not error even without return_best
  expect_no_error(
    result <- comp_to_struc(comps, db = NULL, return_best = FALSE)
  )
  # Result should be a tibble with composition and structure columns
  expect_named(result, c("composition", "structure"))
  # Should have at least one match from default glydb
  expect_gt(nrow(result), 0)
})

test_that("comp_to_struc with return_best=TRUE returns NA for no match", {
  db <- glyrepr::as_glycan_structure(c("Gal(b1-3)GalNAc(a1-"))
  attr(db, "confidence") <- c(1.0)

  comps <- glyrepr::as_glycan_composition(c(Gal = 1, GalNAc = 1))
  comps_no_match <- glyrepr::as_glycan_composition(c(Gal = 5)) # Won't match
  all_comps <- c(comps, comps_no_match)

  result <- comp_to_struc(all_comps, db, return_best = TRUE)

  expect_false(is.data.frame(result))
  expect_equal(length(result), 2)
  expect_equal(as.character(result[1]), "Gal(b1-3)GalNAc(a1-")
  expect_true(is.na(result[2]))
})
