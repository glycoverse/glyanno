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
  expect_equal(
    as.character(res$enhanced),
    c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
  )
})

test_that("enhance_struc enhances partial to intact level with wildcard", {
  db_intact <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
  input_partial <- "Gal(b1-?)GalNAc(a1-"
  res <- enhance_struc(input_partial, to_level = "intact", db = db_intact)
  expect_equal(
    as.character(res$enhanced),
    c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
  )
})

test_that("enhance_struc enhances specific partial to intact level", {
  db_intact <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
  input_partial_specific <- "Gal(b1-3)GalNAc(a1-"
  res <- enhance_struc(
    input_partial_specific,
    to_level = "intact",
    db = db_intact
  )
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

test_that("enhance_struc with return_best=TRUE filters enhanced structures by confidence", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  attr(db, "confidence") <- c(1.0, 2.0)

  struc <- glyrepr::as_glycan_structure("Hex(??-?)HexNAc(??-")
  result <- enhance_struc(
    struc,
    to_level = "intact",
    db = db,
    return_best = TRUE
  ) |>
    dplyr::mutate(enhanced = as.character(enhanced))

  expect_equal(nrow(result), 1)
  expect_equal(result$enhanced, "Gal(b1-4)GalNAc(a1-")
})

test_that("enhance_struc with return_best=TRUE keeps all already-at-level structures", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  attr(db, "confidence") <- c(1.0, 2.0)

  struc <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  result <- enhance_struc(
    struc,
    to_level = "intact",
    db = db,
    return_best = TRUE
  )

  expect_equal(nrow(result), 2)
})

test_that("enhance_struc with return_best=TRUE errors without confidence attr", {
  db <- glyrepr::as_glycan_structure("HexNAc(??-")
  struc <- glyrepr::as_glycan_structure("HexNAc(??-")

  expect_error(
    enhance_struc(struc, to_level = "intact", db = db, return_best = TRUE),
    "must have a .*confidence.* attribute"
  )
})

test_that("enhance_struc treats NA confidence as lowest for enhancement", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  attr(db, "confidence") <- c(NA_real_, 2.0)

  struc <- glyrepr::as_glycan_structure("Hex(??-?)HexNAc(??-")
  result <- enhance_struc(
    struc,
    to_level = "intact",
    db = db,
    return_best = TRUE
  ) |>
    dplyr::mutate(enhanced = as.character(enhanced))

  expect_equal(result$enhanced, "Gal(b1-4)GalNAc(a1-")
})
