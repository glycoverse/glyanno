test_that("enhance_struc enhances basic to topological level", {
  # Hex(??-?)HexNAc should match Gal(??-?)GalNAc
  db_topo <- "Gal(??-?)GalNAc(??-"
  input_basic <- "Hex(??-?)HexNAc(??-"
  res <- enhance_struc(input_basic, db = db_topo)
  expect_equal(as.character(res$enhanced), "Gal(??-?)GalNAc(??-")
})

test_that("enhance_struc enhances topological to intact level", {
  db_intact <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  )
  input_topo <- "Gal(??-?)GalNAc(??-"
  res <- enhance_struc(input_topo, db = db_intact)
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
  res <- enhance_struc(input_partial, db = db_intact)
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
  res <- enhance_struc(input_intact, db = db_intact)
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
  res <- enhance_struc(input_nomatch, db = db_intact)
  expect_equal(res$raw, glyrepr::glycan_structure())
  expect_equal(res$enhanced, glyrepr::glycan_structure())
})

test_that("enhance_struc accepts empty structures", {
  db_intact <- "Gal(b1-3)GalNAc(a1-"
  input_empty <- glyrepr::glycan_structure()

  res <- enhance_struc(input_empty, db = db_intact)

  expect_equal(
    res,
    tibble::tibble(
      raw = glyrepr::glycan_structure(),
      enhanced = glyrepr::glycan_structure()
    )
  )
})

test_that("enhance_struc accepts all-NA structures", {
  db_intact <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  attr(db_intact, "confidence") <- 1
  input_na <- c(glyrepr::glycan_structure(NA), glyrepr::glycan_structure(NA))

  res <- enhance_struc(input_na, db = db_intact)
  best <- enhance_struc(input_na, db = db_intact, return_best = TRUE)

  expect_equal(
    res,
    tibble::tibble(
      raw = input_na,
      enhanced = input_na
    )
  )
  expect_equal(best, input_na)
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
  res <- enhance_struc(input_strucs, db = db_intact)
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
    db = db,
    return_best = TRUE
  )

  expect_equal(length(result), 1)
  expect_equal(as.character(result), "Gal(b1-4)GalNAc(a1-")
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
    db = db,
    return_best = TRUE
  )

  expect_equal(length(result), 2)
})

test_that("enhance_struc with return_best=TRUE errors without confidence attr", {
  db <- glyrepr::as_glycan_structure("HexNAc(??-")
  struc <- glyrepr::as_glycan_structure("HexNAc(??-")

  expect_error(
    enhance_struc(struc, db = db, return_best = TRUE),
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
    db = db,
    return_best = TRUE
  )

  expect_equal(as.character(result), "Gal(b1-4)GalNAc(a1-")
})

test_that("enhance_struc works with db=NULL (regression: confidences undefined)", {
  # This test ensures confidences variable is defined when db=NULL
  struc <- glyrepr::as_glycan_structure("Hex(??-?)HexNAc(??-")
  # Should not error even without return_best
  expect_no_error(
    result <- enhance_struc(
      struc,
      db = NULL,
      return_best = FALSE
    )
  )
  # Result should be a tibble with raw and enhanced columns
  expect_named(result, c("raw", "enhanced"))
  # Should have at least one match from default glydb
  expect_gt(nrow(result), 0)
})

test_that("enhance_struc uses scalar db structure level", {
  db_partial <- c(
    "Gal(??-?)GalNAc(??-",
    "Gal(b1-3)GalNAc(a1-"
  )
  input_basic <- "Hex(??-?)HexNAc(??-"

  res <- enhance_struc(input_basic, db = db_partial)

  expect_equal(
    as.character(res$enhanced),
    c("Gal(??-?)GalNAc(??-", "Gal(b1-3)GalNAc(a1-")
  )
})

test_that("enhance_struc uses scalar input structure level", {
  db_intact <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  attr(db_intact, "confidence") <- 1
  input_partial <- c(
    "Gal(??-?)GalNAc(??-",
    "Gal(b1-4)GalNAc(a1-"
  )

  res <- enhance_struc(input_partial, db = db_intact, return_best = TRUE)

  expect_equal(length(res), 2)
  expect_equal(as.character(res[1]), "Gal(b1-3)GalNAc(a1-")
  expect_true(is.na(res[2]))
})

test_that("enhance_struc with return_best=TRUE returns NA for no match", {
  # db at intact level with two entries
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-"
  ))
  attr(db, "confidence") <- c(1.0, 1.0)

  # First input matches (topological pattern matches intact Gal structures)
  # Second input doesn't match any intact GalNAc structure (Man pattern wrong)
  strucs <- glyrepr::as_glycan_structure(c(
    "Gal(??-?)GalNAc(??-",
    "Man(??-?)Man(??-"
  ))
  result <- enhance_struc(strucs, db = db, return_best = TRUE)

  expect_false(is.data.frame(result))
  expect_equal(length(result), 2)
  expect_equal(as.character(result[1]), "Gal(b1-3)GalNAc(a1-")
  expect_true(is.na(result[2]))
})

test_that("enhance_struc de-novo enhances an N-glycan core", {
  input <- glyrepr::n_glycan_core(
    linkage = FALSE,
    mono_type = "generic"
  )
  expected <- glyrepr::n_glycan_core(
    linkage = FALSE,
    mono_type = "concrete"
  )

  result <- enhance_struc(
    input,
    method = "denovo"
  )

  expect_equal(result$raw, input)
  expect_equal(as.character(result$enhanced), as.character(expected))
})

test_that("enhance_struc de-novo combines all matched branches", {
  input <- paste0(
    "HexNAc(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]Hex(??-?)",
    "HexNAc(??-?)HexNAc(??-"
  )

  result <- enhance_struc(
    input,
    method = "denovo"
  )

  expect_equal(nrow(result), 3)
  expect_equal(length(unique(result$enhanced)), 3)
  expect_equal(
    unique(glyrepr::get_structure_level(result$enhanced)),
    "topological"
  )
})

test_that("enhance_struc de-novo ranks branches and enhances core additions", {
  input <- paste0(
    "HexNAc(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]",
    "[HexNAc(??-?)]Hex(??-?)HexNAc(??-?)[dHex(??-?)]HexNAc(??-"
  )
  expected <- paste0(
    "GlcNAc(??-?)Man(??-?)[GlcNAc(??-?)Man(??-?)]",
    "[GlcNAc(??-?)]Man(??-?)GlcNAc(??-?)[Fuc(??-?)]GlcNAc(??-"
  )

  result <- enhance_struc(
    input,
    method = "denovo",
    return_best = TRUE
  )

  expect_equal(as.character(result), expected)
})

test_that("enhance_struc de-novo handles hybrid N-glycans", {
  input <- paste0(
    "Hex(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]",
    "Hex(??-?)HexNAc(??-?)HexNAc(??-"
  )
  expected <- paste0(
    "Gal(??-?)GlcNAc(??-?)Man(??-?)[Man(??-?)Man(??-?)]",
    "Man(??-?)GlcNAc(??-?)GlcNAc(??-"
  )

  result <- enhance_struc(
    input,
    method = "denovo",
    return_best = TRUE
  )

  expect_equal(as.character(result), expected)
})

test_that("enhance_struc de-novo constrains high-mannose N-glycans", {
  reference <- paste0(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)",
    "[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ) |>
    glyparse::auto_parse()
  input <- glyrepr::reduce_structure_level(reference, "basic")
  expected <- glyrepr::reduce_structure_level(reference, "topological")

  result <- enhance_struc(
    input,
    method = "denovo",
    return_best = TRUE
  )

  expect_equal(as.character(result), as.character(expected))
})

test_that("enhance_struc returns only high-mannose reference subtrees", {
  input <- paste0(
    "Hex(??-?)Hex(??-?)[Hex(??-?)]Hex(??-?)",
    "HexNAc(??-?)HexNAc(??-"
  )
  expected <- paste0(
    "Man(??-?)Man(??-?)[Man(??-?)]Man(??-?)",
    "GlcNAc(??-?)GlcNAc(??-"
  )

  result <- enhance_struc(
    input,
    method = "denovo"
  )

  expect_equal(as.character(result$enhanced), expected)
})
