test_that("comp_to_struc works for generic db and generic comps", {
  db <- glyrepr::as_glycan_structure(c(
    "Hex(??-?)[HexNAc(??-?)]HexNAc(??-",
    "Hex(??-?)HexNAc(??-",
    "HexNAc(??-"
  ))
  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  result <- comp_to_struc(comps, db) |>
    dplyr::mutate(structure = as.character(structure))
  expected <- tibble::tibble(composition = comps, structure = "Hex(??-?)HexNAc(??-")
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

test_that("comp_to_struc works for generic db and concrete comps", {
  # This test has two purposes:
  # 1. Ensure that generic compositions can not match concrete structures in `db`.
  # 2. Ensure that zero-row result has the expected glycan vector types.
  db <- glyrepr::as_glycan_structure("HexNAc(??-")
  comps <- glyrepr::as_glycan_composition("GalNAc(1)")
  result <- comp_to_struc(comps, db)
  expected <- tibble::tibble(
    composition = glyrepr::glycan_composition(),
    structure = glyrepr::glycan_structure()
  )
  expect_equal(result, expected)
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

test_that("comp_to_struc works for mixed db and mixed comps", {
  db <- glyrepr::as_glycan_structure(c(
    "Hex(??-?)HexNAc(??-",  # generic
    "Gal(b1-3)GalNAc(a1-",  # concrete
    "Gal(b1-4)GalNAc(a1-"  # concrete
  ))
  comps <- glyrepr::as_glycan_composition(c(
    "Hex(1)HexNAc(1)",  # generic
    "Gal(1)GalNAc(1)"  # concrete
  ))
  result <- comp_to_struc(comps, db) |>
    dplyr::mutate(structure = as.character(structure))
  expected <- tibble::tibble(
    composition = c(rep(comps[1], 3), rep(comps[2], 2)),
    structure = as.character(c(db, db[2:3]))
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
  comps <- glyrepr::as_glycan_composition(c("Hex(1)HexNAc(1)", "HexNAc(1)", "Hex(1)HexNAc(1)"))
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