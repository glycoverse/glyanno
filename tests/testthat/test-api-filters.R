test_that("structure filters match the corresponding glydb database", {
  withr::local_options(lifecycle_verbosity = "quiet")
  filters <- list(
    glycan_type = "O-GalNAc",
    species = "Homo sapiens",
    structure_level = "topological",
    mono_type = "concrete",
    mono_range = list(Gal = c(1L, 1L), GalNAc = c(1L, 1L))
  )
  db <- do.call(glydb::glydb_structures, filters)
  comp <- "Gal(1)GalNAc(1)"
  struc <- "Gal(??-?)GalNAc(??-"

  expected_comp <- suppressWarnings(comp_to_struc(comp, db = db))
  actual_comp <- suppressWarnings(do.call(
    comp_to_struc,
    c(list(comps = comp), filters)
  ))
  expect_equal(actual_comp, expected_comp)

  expected_struc <- suppressWarnings(enhance_struc(struc, db = db))
  actual_struc <- suppressWarnings(do.call(
    enhance_struc,
    c(list(strucs = struc), filters)
  ))
  expect_identical(
    as.character(actual_struc$raw),
    as.character(expected_struc$raw)
  )
  expect_identical(
    as.character(actual_struc$enhanced),
    as.character(expected_struc$enhanced)
  )
  expect_equal(actual_struc$confidence, expected_struc$confidence)
})

test_that("composition filters match the corresponding glydb database", {
  withr::local_options(lifecycle_verbosity = "quiet")
  filters <- list(
    glycan_type = "O-GalNAc",
    species = "Homo sapiens",
    mono_type = "concrete",
    mono_range = list(Gal = c(1L, 1L), GalNAc = c(1L, 1L))
  )
  db <- do.call(glydb::glydb_compositions, filters)
  expected_mz <- mz_to_comp(406.1325, adduct = "Na+", db = db)
  actual_mz <- do.call(
    mz_to_comp,
    c(list(mz = 406.1325, adduct = "Na+"), filters)
  )
  expect_equal(actual_mz, expected_mz)

  expected_enhanced <- enhance_comp("Hex(1)HexNAc(1)", db = db)
  actual_enhanced <- do.call(
    enhance_comp,
    c(list(comps = "Hex(1)HexNAc(1)"), filters)
  )
  expect_equal(actual_enhanced, expected_enhanced)
})

test_that("species filter accepts the capitalization in the new API example", {
  lower <- suppressWarnings(comp_to_struc(
    "H5N2",
    glycan_type = "N",
    species = "Homo sapiens",
    return_best = TRUE
  ))
  title <- suppressWarnings(comp_to_struc(
    "H5N2",
    glycan_type = "N",
    species = "Homo Sapiens",
    return_best = TRUE
  ))
  expect_equal(title, lower)
  expect_length(title, 1L)
  expect_equal(is.na(title), FALSE)
})

test_that("deprecated db remains usable and rejects mixed filtering", {
  db <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  lifecycle::expect_deprecated(comp_to_struc("Gal(1)GalNAc(1)", db = db))
  lifecycle::expect_deprecated(enhance_struc("Gal(??-?)GalNAc(??-", db = db))
  comp_db <- glyrepr::as_glycan_composition("Gal(1)GalNAc(1)")
  lifecycle::expect_deprecated(enhance_comp("Hex(1)HexNAc(1)", db = comp_db))
  lifecycle::expect_deprecated(mz_to_comp(
    406.1325,
    adduct = "Na+",
    db = comp_db
  ))
  expect_snapshot(error = TRUE, {
    suppressWarnings(comp_to_struc(
      "Gal(1)GalNAc(1)",
      db = db,
      glycan_type = "O-GalNAc"
    ))
  })
})
