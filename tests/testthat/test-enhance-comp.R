test_that("enhance_comp works for a basic case", {
  comps <- glyrepr::glycan_composition(c(Hex = 1, HexNAc = 1))
  db <- glyrepr::glycan_composition(c(Gal = 1, GalNAc = 1), c(Glc = 1, GalNAc = 1))
  result <- enhance_comp(comps, db)
  expected <- tibble::tibble(raw = comps, enhanced = db)
  expect_equal(result, expected)
})

test_that("enhance_comp accepts character input", {
  comps <- "Hex(1)HexNAc(1)"
  db <- c("Gal(1)GalNAc(1)", "Glc(1)GalNAc(1)")
  result <- enhance_comp(comps, db)
  expected <- tibble::tibble(
    raw = glyrepr::glycan_composition(c(Hex = 1, HexNAc = 1)),
    enhanced = glyrepr::glycan_composition(c(Gal = 1, GalNAc = 1), c(Glc = 1, GalNAc = 1))
  )
  expect_equal(result, expected)
})

test_that("enhance_comp works with empty input", {
  comps <- glyrepr::glycan_composition()
  db <- glyrepr::glycan_composition(c(Gal = 1, GalNAc = 1), c(Glc = 1, GalNAc = 1))
  result <- enhance_comp(comps, db)
  expected <- tibble::tibble(raw = glyrepr::glycan_composition(), enhanced = glyrepr::glycan_composition())
  expect_equal(result, expected)
})

test_that("enhance_comp works with generic and compositions", {
  comps <- "Hex(1)HexNAc(1)"
  db <- c("Glc(1)GalNAc(1)", "Glc(1)GlcNAc(1)")
  result <- enhance_comp(comps, db)
  expected <- tibble::tibble(
    raw = glyrepr::as_glycan_composition(c("Hex(1)HexNAc(1)", "Hex(1)HexNAc(1)")),
    enhanced = glyrepr::as_glycan_composition(c("Glc(1)GalNAc(1)", "Glc(1)GlcNAc(1)"))
  )
  expect_equal(result, expected)
})

test_that("enhance_comp works with concrete compositions", {
  comps <- "Gal(1)GalNAc(1)"
  db <- c("Glc(1)GalNAc(1)", "Glc(1)GlcNAc(1)")
  result <- enhance_comp(comps, db)
  expected <- tibble::tibble(
    raw = glyrepr::as_glycan_composition("Gal(1)GalNAc(1)"),
    enhanced = glyrepr::as_glycan_composition("Gal(1)GalNAc(1)")
  )
  expect_equal(result, expected)
})

test_that("enhance_comp issues a warning when some compositions in db are generic", {
  comps <- "Hex(1)HexNAc(1)"
  db <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  expect_warning(enhance_comp(comps, db), "Some compositions in `db` are generic, which will be dropped.")
})