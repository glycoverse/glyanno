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
    structure = "Hex(??-?)HexNAc(??-",
    confidence = NA_real_
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
    structure = c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-"),
    confidence = c(NA_real_, NA_real_)
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc leaves incompatible compositions unmatched", {
  db <- glyrepr::as_glycan_structure("HexNAc(??-")
  comps <- glyrepr::as_glycan_composition("GalNAc(1)")

  expanded <- comp_to_struc(comps, db)
  attr(db, "confidence") <- 1
  best <- comp_to_struc(comps, db, return_best = TRUE)

  expect_equal(
    expanded,
    tibble::tibble(
      composition = glyrepr::glycan_composition(),
      structure = glyrepr::glycan_structure(),
      confidence = numeric()
    )
  )
  expect_identical(is.na(best), TRUE)
})

test_that("comp_to_struc matches mixed compositions residue by residue", {
  db <- glyrepr::as_glycan_structure(c(
    "Man(??-?)Man(??-",
    "Gal(??-?)Man(??-",
    "Hex(??-?)Man(??-",
    "Hex(??-?)Hex(??-"
  ))
  attr(db, "confidence") <- 1:4
  comps <- glyrepr::as_glycan_composition(c(
    "Man(1)Hex(1)",
    "Man(2)",
    "Hex(2)",
    NA
  ))

  expanded <- comp_to_struc(comps, db)
  best <- comp_to_struc(comps, db, return_best = TRUE)

  expect_equal(
    as.character(expanded$structure),
    c(
      "Man(??-?)Man(??-",
      "Gal(??-?)Man(??-",
      "Hex(??-?)Man(??-",
      "Man(??-?)Man(??-",
      "Man(??-?)Man(??-",
      "Gal(??-?)Man(??-",
      "Hex(??-?)Man(??-",
      "Hex(??-?)Hex(??-"
    )
  )
  expect_equal(
    as.character(best),
    c("Hex(??-?)Man(??-", "Man(??-?)Man(??-", "Hex(??-?)Hex(??-", NA)
  )
})

test_that("comp_to_struc excludes floating database structures", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(??-?)Man(??-?)[Man(??-?)]GlcNAc(??-",
    "{Gal(??-?)|2,3}Man(??-?)[Man(??-?)]GlcNAc(??-"
  ))
  attr(db, "confidence") <- c(1, 2)
  comps <- glyrepr::as_glycan_composition("Hex(3)HexNAc(1)")

  expect_snapshot(
    result <- comp_to_struc(comps, db, return_best = TRUE)
  )

  expect_equal(
    as.character(result),
    "Gal(??-?)Man(??-?)[Man(??-?)]GlcNAc(??-"
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
    structure = c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-"),
    confidence = c(NA_real_, NA_real_)
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
    structure = c("Hex(??-?)HexNAc(??-", "HexNAc(??-"),
    confidence = c(NA_real_, NA_real_)
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
    structure = c("Hex(??-?)HexNAc(??-", "HexNAc(??-"),
    confidence = c(NA_real_, NA_real_)
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
  attr(db, "confidence") <- c(1, 2)
  result <- comp_to_struc(comps, db) |>
    dplyr::mutate(structure = as.character(structure))
  best <- comp_to_struc(comps, db, return_best = TRUE)
  expected <- tibble::tibble(
    composition = comps,
    structure = c("Hex(??-?)HexNAc(??-", "HexNAc(??-", "Hex(??-?)HexNAc(??-"),
    confidence = c(1, 2, 1)
  )
  expect_equal(result, expected)
  expect_equal(as.character(best), expected$structure)
})

test_that("comp_to_struc matches each unique composition once", {
  matched <- character()
  local_mocked_bindings(
    .composition_match_ids = function(pattern, index) {
      matched <<- c(matched, as.character(pattern))
      integer()
    },
    .package = "glyanno"
  )
  db <- glyrepr::as_glycan_structure("GalNAc(??-")
  comps <- glyrepr::as_glycan_composition(c(
    "HexNAc(1)",
    "Hex(1)",
    "HexNAc(1)"
  ))

  comp_to_struc(comps, db)

  expect_equal(matched, c("HexNAc(1)", "Hex(1)"))
})

test_that("comp_to_struc reuses bundled metadata for the default database", {
  old_db_cache <- as.list(.default_struc_db_cache, all.names = TRUE)
  old_comp_cache <- as.list(.comp_to_struc_cache, all.names = TRUE)
  on.exit(
    {
      rm(list = ls(.default_struc_db_cache), envir = .default_struc_db_cache)
      rm(list = ls(.comp_to_struc_cache), envir = .comp_to_struc_cache)
      list2env(old_db_cache, envir = .default_struc_db_cache)
      list2env(old_comp_cache, envir = .comp_to_struc_cache)
    },
    add = TRUE
  )
  rm(list = ls(.default_struc_db_cache), envir = .default_struc_db_cache)
  rm(list = ls(.comp_to_struc_cache), envir = .comp_to_struc_cache)

  metadata <- default_comp_to_struc_metadata
  eligible <- which(!metadata$floating)
  distinct <- !duplicated(metadata$generic_keys[eligible])
  ids <- eligible[distinct][1:2]
  db <- glyrepr::as_glycan_structure(metadata$structure_keys[ids])
  attr(db, "confidence") <- c(1, 2)
  inputs <- c(
    glyrepr::convert_to_generic(metadata$composition[ids[1]]),
    metadata$composition[ids[2]]
  )
  local_mocked_bindings(
    glydb_structures = function(...) db,
    .package = "glydb"
  )
  local_mocked_bindings(
    as_glycan_composition = function(...) {
      stop("default structures were converted again")
    },
    .package = "glyrepr"
  )

  result <- comp_to_struc(inputs, return_best = TRUE)

  expect_equal(as.character(result), as.character(db))
})

test_that("default database metadata falls back for unknown structures", {
  db <- glyrepr::as_glycan_structure("Hex(??-?)HexNAc(??-")

  result <- .prepare_default_struc_db(db, "db")

  expect_equal(as.character(result$composition), "Hex(1)HexNAc(1)")
  expect_equal(result$generic_keys, "Hex(1)HexNAc(1)")
})

test_that("comp_to_struc accepts empty compositions", {
  db <- glyrepr::as_glycan_structure(c("Hex(??-?)HexNAc(??-", "HexNAc(??-"))
  comps <- glyrepr::glycan_composition()
  result <- comp_to_struc(comps, db)
  expected <- tibble::tibble(
    composition = glyrepr::glycan_composition(),
    structure = glyrepr::glycan_structure(),
    confidence = numeric()
  )
  expect_equal(result, expected)
})

test_that("comp_to_struc rejects empty db", {
  db <- glyrepr::glycan_structure()
  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")
  expect_error(comp_to_struc(comps, db))
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
    result <- suppressWarnings(
      comp_to_struc(comps, db = NULL, return_best = FALSE)
    )
  )
  expect_named(result, c("composition", "structure", "confidence"))
  expect_false(all(is.na(result$confidence)))
  # Should have at least one match from default glydb
  expect_gt(nrow(result), 0)
})

test_that("comp_to_struc reuses the prepared default database", {
  comps <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)")

  first <- comp_to_struc(comps, db = NULL, return_best = TRUE)
  second <- comp_to_struc(comps, db = NULL, return_best = TRUE)

  expect_equal(second, first)
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
