test_that("fill_anomer_pos fills missing anomer positions", {
  strucs <- glyrepr::as_glycan_structure(c(
    "Gal(??-?)GalNAc(??-",
    "Neu5Ac(??-?)Gal(??-?)GalNAc(??-"
  ))

  result <- suppressMessages(fill_anomer_pos(strucs))

  expect_s3_class(result, "glyrepr_structure")
  expect_equal(
    as.character(result),
    c(
      "Gal(?1-?)GalNAc(?1-",
      "Neu5Ac(?2-?)Gal(?1-?)GalNAc(?1-"
    )
  )
})

test_that("fill_anomer_pos preserves existing anomer annotations", {
  strucs <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Neu5Ac(a?-3)Gal(?1-"
  ))

  result <- suppressMessages(fill_anomer_pos(strucs))

  expect_equal(
    as.character(result),
    c(
      "Gal(b1-3)GalNAc(a1-",
      "Neu5Ac(a2-3)Gal(?1-"
    )
  )
})

test_that("fill_anomer_pos preserves NA values and names", {
  strucs <- c(
    missing = glyrepr::glycan_structure(NA),
    present = glyrepr::as_glycan_structure("Gal(??-?)GalNAc(??-")
  )

  result <- suppressMessages(fill_anomer_pos(strucs))

  expect_equal(names(result), names(strucs))
  expect_true(is.na(result[[1]]))
  expect_equal(as.character(result[[2]]), "Gal(?1-?)GalNAc(?1-")
})

test_that("fill_anomer_pos throws error for non-concrete monosaccharides", {
  strucs <- glyrepr::n_glycan_core(linkage = FALSE, mono_type = "generic")

  expect_error(
    fill_anomer_pos(strucs),
    regexp = "must have concrete monosaccharides"
  )
})
