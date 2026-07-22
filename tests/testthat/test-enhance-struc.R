test_that("enhance_struc enhances basic to topological level", {
  # Hex(??-?)HexNAc should match Gal(??-?)GalNAc
  db_topo <- "Gal(??-?)GalNAc(??-"
  input_basic <- "Hex(??-?)HexNAc(??-"
  res <- enhance_struc(input_basic, db = db_topo)
  expect_equal(as.character(res$enhanced), "Gal(??-?)GalNAc(??-")
})

test_that("structure enhancement functions expose separate contracts", {
  expect_identical(
    names(formals(enhance_struc)),
    c("strucs", "db", "return_best")
  )
  expect_identical(
    names(formals(enhance_struc_denovo)),
    c("strucs", "fallback_db")
  )
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

test_that("enhance_struc restores repeated and missing inputs", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-",
    "GalNAc(a1-"
  ))
  attr(db, "confidence") <- c(1, 2, 3)
  strucs <- glyrepr::as_glycan_structure(c(
    "Hex(??-?)HexNAc(??-",
    NA,
    "HexNAc(??-",
    "Hex(??-?)HexNAc(??-"
  ))

  best <- enhance_struc(strucs, db = db, return_best = TRUE)
  expanded <- enhance_struc(strucs, db = db, return_best = FALSE)

  expect_equal(
    as.character(best),
    c("Gal(b1-4)GalNAc(a1-", NA, "GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
  )
  expect_equal(
    as.character(expanded$raw),
    c(
      "Hex(??-?)HexNAc(??-",
      "Hex(??-?)HexNAc(??-",
      "HexNAc(??-",
      "Hex(??-?)HexNAc(??-",
      "Hex(??-?)HexNAc(??-"
    )
  )
  expect_equal(
    as.character(expanded$enhanced),
    c(
      "Gal(b1-3)GalNAc(a1-",
      "Gal(b1-4)GalNAc(a1-",
      "GalNAc(a1-",
      "Gal(b1-3)GalNAc(a1-",
      "Gal(b1-4)GalNAc(a1-"
    )
  )
})

test_that("enhance_struc matches each unique non-missing input once", {
  matched_patterns <- NULL
  local_mocked_bindings(
    have_motifs = function(glycans, motifs, ...) {
      matched_patterns <<- motifs
      matrix(TRUE, nrow = length(glycans), ncol = length(motifs))
    },
    .package = "glymotif"
  )
  db <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  attr(db, "confidence") <- 1
  strucs <- glyrepr::as_glycan_structure(c(
    rep("Hex(??-?)HexNAc(??-", 100),
    "HexNAc(??-",
    NA
  ))

  enhance_struc(strucs, db = db, return_best = TRUE)

  expect_equal(
    as.character(matched_patterns),
    c("Hex(??-?)HexNAc(??-", "HexNAc(??-")
  )
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

test_that("enhance_struc_denovo enhances an N-glycan core", {
  input <- glyrepr::n_glycan_core(
    linkage = FALSE,
    mono_type = "generic"
  )
  expected <- glyrepr::n_glycan_core(
    linkage = FALSE,
    mono_type = "concrete"
  )

  result <- enhance_struc_denovo(
    input
  )

  expect_equal(as.character(result), as.character(expected))
})

test_that("enhance_struc_denovo keeps one branch per basic pattern", {
  input <- paste0(
    "HexNAc(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]Hex(??-?)",
    "HexNAc(??-?)HexNAc(??-"
  )

  result <- enhance_struc_denovo(
    input
  )

  expect_length(result, 1)
  expect_equal(
    glyrepr::get_structure_level(result),
    "topological"
  )
})

test_that("de-novo branch data contain all sialic acid variants", {
  basic_branches <- glyrepr::reduce_structure_level(
    topological_branches,
    "basic"
  )
  direct_keys <- vapply(
    as.list(basic_branches),
    glyrepr::graph_to_iupac,
    character(1)
  )

  branch_keys <- as.character(topological_branches)
  expect_identical(anyDuplicated(branch_keys), 0L)
  expect_identical(unname(direct_keys), unname(as.character(basic_branches)))
  expect_length(topological_branch_index, length(unique(direct_keys)))
  expect_setequal(
    unname(unlist(topological_branch_index)),
    seq_along(topological_branches)
  )
  expect_equal(
    unname(lengths(topological_branch_index)),
    rep(1L, length(topological_branch_index))
  )
  expected_templates <- purrr::map(
    seq_along(topological_branches),
    \(id) .topological_subtree_template(topological_branches[id])
  )
  names(expected_templates) <- as.character(seq_along(expected_templates))
  expect_equal(topological_branch_templates, expected_templates)
  expect_null(attr(topological_branches, "confidence"))

  for (graph in as.list(topological_branches)) {
    monos <- igraph::vertex_attr(graph, "mono")
    for (node in which(monos %in% c("Neu5Ac", "Neu5Gc"))) {
      replacement <- if (monos[[node]] == "Neu5Ac") "Neu5Gc" else "Neu5Ac"
      variant <- igraph::set_vertex_attr(
        graph,
        "mono",
        index = node,
        value = replacement
      )
      expect_contains(
        branch_keys,
        as.character(glyrepr::as_glycan_structure(variant))
      )
    }
  }

  expect_identical(
    exists(
      "intact_branches",
      envir = asNamespace("glyanno"),
      inherits = FALSE
    ),
    FALSE
  )
})

test_that("enhance_struc_denovo ranks branches and enhances core additions", {
  input <- paste0(
    "HexNAc(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]",
    "[HexNAc(??-?)]Hex(??-?)HexNAc(??-?)[dHex(??-?)]HexNAc(??-"
  )
  expected <- paste0(
    "GlcNAc(??-?)Man(??-?)[GlcNAc(??-?)Man(??-?)]",
    "[GlcNAc(??-?)]Man(??-?)GlcNAc(??-?)[Fuc(??-?)]GlcNAc(??-"
  )

  result <- enhance_struc_denovo(
    input
  )

  expect_equal(as.character(result), expected)
})

test_that("de-novo table materialization follows the validated graph pipeline", {
  expected <- glyrepr::as_glycan_structure(paste0(
    "GlcNAc(??-?)Man(??-?)[GlcNAc(??-?)Man(??-?)]",
    "Man(??-?)GlcNAc(??-?)GlcNAc(??-"
  ))
  candidate <- list(
    nodes = glyrepr::structure_nodes(expected),
    edges = glyrepr::structure_edges(expected)
  )

  result <- .topological_n_glycan_from_tables(list(candidate, candidate))

  expect_s3_class(result, "glyrepr_structure")
  expect_identical(as.character(result), rep(as.character(expected), 2))
  expect_length(attr(result, "graphs"), 1)
})

test_that("enhance_struc_denovo handles hybrid N-glycans", {
  input <- paste0(
    "Hex(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]",
    "Hex(??-?)HexNAc(??-?)HexNAc(??-"
  )
  expected <- paste0(
    "Gal(??-?)GlcNAc(??-?)Man(??-?)[Man(??-?)Man(??-?)]",
    "Man(??-?)GlcNAc(??-?)GlcNAc(??-"
  )

  result <- enhance_struc_denovo(
    input
  )

  expect_equal(as.character(result), expected)
})

test_that("enhance_struc_denovo constrains high-mannose N-glycans", {
  reference <- paste0(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)",
    "[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ) |>
    glyparse::auto_parse()
  input <- glyrepr::reduce_structure_level(reference, "basic")
  expected <- glyrepr::reduce_structure_level(reference, "topological")

  result <- enhance_struc_denovo(
    input
  )

  expect_equal(as.character(result), as.character(expected))
})

test_that("enhance_struc_denovo returns high-mannose reference subtrees", {
  input <- paste0(
    "Hex(??-?)Hex(??-?)[Hex(??-?)]Hex(??-?)",
    "HexNAc(??-?)HexNAc(??-"
  )
  expected <- paste0(
    "Man(??-?)Man(??-?)[Man(??-?)]Man(??-?)",
    "GlcNAc(??-?)GlcNAc(??-"
  )

  result <- enhance_struc_denovo(
    input
  )

  expect_equal(as.character(result), expected)
})

test_that("enhance_struc_denovo preserves high-mannose core additions", {
  input <- paste0(
    "Hex(??-?)Hex(??-?)[Hex(??-?)][HexNAc(??-?)]Hex(??-?)",
    "HexNAc(??-?)[dHex(??-?)]HexNAc(??-"
  )
  expected <- paste0(
    "Man(??-?)Man(??-?)[Man(??-?)][GlcNAc(??-?)]Man(??-?)",
    "GlcNAc(??-?)[Fuc(??-?)]GlcNAc(??-"
  )

  result <- enhance_struc_denovo(
    input
  )

  expect_equal(as.character(result), expected)
})

test_that("enhance_struc_denovo falls back to database matching", {
  db <- glyrepr::as_glycan_structure("Gal(??-?)GlcNAc(??-")
  attr(db, "confidence") <- 1
  input <- c(
    "Hex(??-?)HexNAc(??-",
    "Hex(??-?)Hex(??-"
  )

  result <- enhance_struc_denovo(
    input,
    fallback_db = db
  )

  expect_equal(as.character(result), c(as.character(db), NA_character_))
})

test_that("enhance_struc_denovo defaults to a topological fallback", {
  result <- enhance_struc_denovo(
    "Hex(??-?)HexNAc(??-"
  )

  expect_length(result, 1)
  expect_equal(
    glyrepr::get_structure_level(result),
    "topological"
  )
})

test_that("enhance_struc_denovo rejects a basic fallback database", {
  input <- glyrepr::n_glycan_core(
    linkage = FALSE,
    mono_type = "generic"
  )

  expect_snapshot(
    error = TRUE,
    enhance_struc_denovo(
      input,
      fallback_db = input
    )
  )
})

test_that("enhance_struc_denovo reduces fallback databases", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-?)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-",
    "Man(b1-4)GlcNAc(b1-"
  ))
  attr(db, "confidence") <- c(1, 2, 3)
  input <- "Hex(??-?)HexNAc(??-"

  result <- enhance_struc_denovo(
    input,
    fallback_db = db
  )

  expect_equal(glyrepr::get_structure_level(result), "topological")
  expect_equal(as.character(result), "Man(??-?)GlcNAc(??-")
})

test_that("de-novo fallback preparation preserves all-missing confidence", {
  db <- glyrepr::as_glycan_structure(c(
    "Gal(b1-?)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-"
  ))
  attr(db, "confidence") <- c(NA_real_, NA_real_)

  result <- .prepare_denovo_struc_db(db)

  expect_length(result, 1)
  expect_identical(attr(result, "confidence"), NA_real_)
})

test_that("enhance_struc_denovo falls back after deduction errors", {
  input <- paste0(
    "dHex(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]Hex(??-?)",
    "HexNAc(??-?)HexNAc(??-"
  )
  db <- glyrepr::as_glycan_structure(paste0(
    "Fuc(??-?)Man(??-?)[GlcNAc(??-?)Man(??-?)]Man(??-?)",
    "GlcNAc(??-?)GlcNAc(??-"
  ))
  attr(db, "confidence") <- 1

  result <- enhance_struc_denovo(
    input,
    fallback_db = db
  )

  expect_equal(as.character(result), as.character(db))
})

test_that("enhance_struc_denovo preserves repeated inputs", {
  input <- paste0(
    "HexNAc(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]Hex(??-?)",
    "HexNAc(??-?)HexNAc(??-"
  )

  single <- enhance_struc_denovo(
    input
  )
  repeated <- enhance_struc_denovo(
    rep(input, 2)
  )

  expect_equal(as.character(repeated), rep(as.character(single), 2))
})

test_that("enhance_struc_denovo batches distinct core matches", {
  inputs <- c(
    paste0(
      "HexNAc(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]Hex(??-?)",
      "HexNAc(??-?)HexNAc(??-"
    ),
    paste0(
      "Hex(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]",
      "Hex(??-?)HexNAc(??-?)HexNAc(??-"
    ),
    "Hex(??-?)HexNAc(??-"
  )
  db <- glyrepr::as_glycan_structure("Gal(??-?)GlcNAc(??-")
  attr(db, "confidence") <- 1

  batched <- enhance_struc_denovo(
    inputs,
    fallback_db = db
  )
  individual <- vapply(
    inputs,
    function(input) {
      as.character(enhance_struc_denovo(
        input,
        fallback_db = db
      ))
    },
    character(1)
  )

  expect_identical(as.character(batched), unname(individual))
})

test_that("enhance_struc_denovo assembles repeated and missing results", {
  input <- glyrepr::n_glycan_core(
    linkage = FALSE,
    mono_type = "generic"
  )
  expected <- glyrepr::n_glycan_core(
    linkage = FALSE,
    mono_type = "concrete"
  )

  result <- enhance_struc_denovo(
    c(input, NA, input)
  )

  expect_s3_class(result, "glyrepr_structure")
  expect_null(names(result))
  expect_equal(
    as.character(result),
    c(as.character(expected), NA_character_, as.character(expected))
  )
})

test_that("enhance_struc_denovo batches repeated database fallbacks", {
  input <- rep("Hex(??-?)HexNAc(??-", 2)
  db <- glyrepr::as_glycan_structure("Gal(??-?)GlcNAc(??-")
  attr(db, "confidence") <- 1

  result <- enhance_struc_denovo(
    input,
    fallback_db = db
  )

  expect_equal(as.character(result), rep(as.character(db), 2))
})
