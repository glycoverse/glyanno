test_that("enhance_struc enhances basic to topological level", {
  # Hex(??-?)HexNAc should match Gal(??-?)GalNAc
  db_topo <- "Gal(??-?)GalNAc(??-"
  input_basic <- "Hex(??-?)HexNAc(??-"
  res <- enhance_struc(input_basic, to_level = "topological", db = db_topo)
  expect_equal(as.character(res$enhanced), "Gal(??-?)GalNAc(??-")
})

test_that("enhance_struc enhances topological to intact level", {
  db_intact <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
  input_topo <- "Gal(??-?)GalNAc(??-"
  res <- enhance_struc(input_topo, to_level = "intact", db = db_intact)
  # Should match both
  expect_equal(as.character(res$enhanced), c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-"))
})

test_that("enhance_struc enhances partial to intact level with wildcard", {
  db_intact <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
  input_partial <- "Gal(b1-?)GalNAc(a1-"
  res <- enhance_struc(input_partial, to_level = "intact", db = db_intact)
  expect_equal(as.character(res$enhanced), c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-"))
})

test_that("enhance_struc enhances specific partial to intact level", {
  db_intact <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
  input_partial_specific <- "Gal(b1-3)GalNAc(a1-"
  res <- enhance_struc(input_partial_specific, to_level = "intact", db = db_intact)
  expect_equal(as.character(res$enhanced), "Gal(b1-3)GalNAc(a1-")
})

test_that("enhance_struc returns unchanged when no enhancement needed", {
  db_intact <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
  input_intact <- "Gal(b1-3)GalNAc(a1-"
  res <- enhance_struc(input_intact, to_level = "intact", db = db_intact)
  expect_equal(as.character(res$enhanced), input_intact)
  expect_equal(as.character(res$raw), input_intact)
})

test_that("enhance_struc returns empty result when no match found", {
  db_intact <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
  # Use a lower level structure that is not in DB
  input_nomatch <- "Man(??-?)Man(??-"
  res <- enhance_struc(input_nomatch, to_level = "intact", db = db_intact)
  expect_equal(res$raw, glyrepr::glycan_structure())
  expect_equal(res$enhanced, glyrepr::glycan_structure())
})

test_that("enhance_struc works with multiple glycans", {
  db_intact <- c(
    "GalNAc(a1-",
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
  input_strucs <- c(
    "Gal(??-?)GalNAc(??-",
    "GalNAc(??-"
  )
  res <- enhance_struc(input_strucs, to_level = "intact", db = db_intact)
  expected <- tibble::tibble(
    raw = glyrepr::as_glycan_structure(c(
      "Gal(??-?)GalNAc(??-",
      "Gal(??-?)GalNAc(??-",
      "GalNAc(??-"
    )),
    enhanced = glyrepr::as_glycan_structure(c(
      "Gal(b1-3)GalNAc(a1-",
      "Gal(b1-4)GalNAc(a1-",
      "GalNAc(a1-"
    ))
  )
  expect_equal(as.character(res$raw), as.character(expected$raw))
  expect_equal(as.character(res$enhanced), as.character(expected$enhanced))
})