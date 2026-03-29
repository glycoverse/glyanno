test_that("enhance_comp works for a basic case", {
  comps <- glyrepr::glycan_composition(c(Hex = 1, HexNAc = 1))
  db <- glyrepr::glycan_composition(
    c(Gal = 1, GalNAc = 1),
    c(Glc = 1, GalNAc = 1)
  )
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
    enhanced = glyrepr::glycan_composition(
      c(Gal = 1, GalNAc = 1),
      c(Glc = 1, GalNAc = 1)
    )
  )
  expect_equal(result, expected)
})

test_that("enhance_comp works with empty input", {
  comps <- glyrepr::glycan_composition()
  db <- glyrepr::glycan_composition(
    c(Gal = 1, GalNAc = 1),
    c(Glc = 1, GalNAc = 1)
  )
  result <- enhance_comp(comps, db)
  expected <- tibble::tibble(
    raw = glyrepr::glycan_composition(),
    enhanced = glyrepr::glycan_composition()
  )
  expect_equal(result, expected)
})

test_that("enhance_comp works with generic and compositions", {
  comps <- "Hex(1)HexNAc(1)"
  db <- c("Glc(1)GalNAc(1)", "Glc(1)GlcNAc(1)")
  result <- enhance_comp(comps, db)
  expected <- tibble::tibble(
    raw = glyrepr::as_glycan_composition(c(
      "Hex(1)HexNAc(1)",
      "Hex(1)HexNAc(1)"
    )),
    enhanced = glyrepr::as_glycan_composition(c(
      "Glc(1)GalNAc(1)",
      "Glc(1)GlcNAc(1)"
    ))
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

test_that("enhance_comp errors when db has generic compositions", {
  comps <- "Hex(1)HexNAc(1)"
  db <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)") # generic composition
  expect_error(
    enhance_comp(comps, db),
    "All compositions in `db` must be concrete"
  )
})

test_that("enhance_comp with return_best=TRUE keeps highest confidence match", {
  db <- glyrepr::glycan_composition(
    c(Glc = 2),
    c(Glc = 1, Gal = 1)
  )
  attr(db, "confidence") <- c(1.0, 2.0)

  comps <- glyrepr::as_glycan_composition("Hex(2)")
  result <- enhance_comp(comps, db, return_best = TRUE)

  expect_equal(length(result), 1)
  expect_equal(as.character(result), "Glc(1)Gal(1)")
})

test_that("enhance_comp with return_best=TRUE errors without confidence attr", {
  db <- glyrepr::glycan_composition(c(Glc = 1))
  comps <- glyrepr::as_glycan_composition("Hex(1)")

  expect_error(
    enhance_comp(comps, db, return_best = TRUE),
    "must have a .*confidence.* attribute"
  )
})

test_that("enhance_comp treats NA confidence as lowest", {
  db <- glyrepr::glycan_composition(
    c(Glc = 2),
    c(Glc = 1, Gal = 1)
  )
  attr(db, "confidence") <- c(NA_real_, 2.0)

  comps <- glyrepr::as_glycan_composition("Hex(2)")
  result <- enhance_comp(comps, db, return_best = TRUE)

  expect_equal(as.character(result), "Glc(1)Gal(1)")
})

test_that("enhance_comp works with db=NULL (regression: confidences undefined)", {
  # This test ensures confidences variable is defined when db=NULL
  comps <- glyrepr::as_glycan_composition("Hex(1)")
  # Should not error even without return_best
  expect_no_error(
    result <- enhance_comp(comps, db = NULL, return_best = FALSE)
  )
  # Result should be a tibble with raw and enhanced columns
  expect_named(result, c("raw", "enhanced"))
})

test_that("enhance_comp with return_best=TRUE returns vector with NA for no match", {
  db <- glyrepr::glycan_composition(
    c(Glc = 2)
  )
  attr(db, "confidence") <- c(1.0)

  comps <- glyrepr::as_glycan_composition(c("Hex(2)", "Hex(3)"))
  result <- enhance_comp(comps, db, return_best = TRUE)

  # Should return a vector (glyrepr::glycan_composition()), not tibble
  expect_false(is.data.frame(result))
  # Length should match input
  expect_equal(length(result), length(comps))
  # First should match, second should be NA
  expect_equal(as.character(result[1]), "Glc(2)")
  expect_true(is.na(result[2]))
})

test_that("enhance_comp with return_best=TRUE and concrete comp returns vector", {
  db <- glyrepr::glycan_composition(
    c(Glc = 1, Gal = 1)
  )
  attr(db, "confidence") <- c(1.0)

  comps <- glyrepr::as_glycan_composition("Glc(1)Gal(1)")
  result <- enhance_comp(comps, db, return_best = TRUE)

  # Should return a vector
  expect_false(is.data.frame(result))
  expect_equal(length(result), 1)
  expect_equal(as.character(result), "Glc(1)Gal(1)")
})
