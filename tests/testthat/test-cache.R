local_cache_fixture <- function(env = parent.frame()) {
  withr::local_options(
    glyanno.cache_dir = tempfile("glyanno-cache-"),
    .local_envir = env
  )
  withr::defer(
    unlink(getOption("glyanno.cache_dir"), recursive = TRUE),
    envir = env
  )
  local_mocked_bindings(
    .annotation_cache = new.env(parent = emptyenv()),
    .env = env
  )
  comp <- glyrepr::as_glycan_composition(c("Gal(1)GalNAc(1)", "Glc(2)"))
  struc <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Glc(a1-4)Glc(b1-"
  ))
  attr(comp, "confidence") <- c(1, 2)
  attr(struc, "confidence") <- c(1, 2)
  local_mocked_bindings(
    glydb_compositions = function(...) comp,
    glydb_structures = function(...) struc,
    .package = "glydb",
    .env = env
  )
  local_mocked_bindings(
    .build_accession_lookup = function() c(test = "G00001"),
    .env = env
  )
}

annotation_semantics <- function(x) {
  if (glyrepr::is_glycan_structure(x)) {
    graphs <- as.list(x[!is.na(x) & !duplicated(as.character(x))])
    return(list(
      values = as.character(x),
      names = names(x),
      confidence = attr(x, "confidence"),
      graphs = lapply(graphs, function(g) {
        list(
          vertices = igraph::vertex_attr(g),
          edges = igraph::as_edgelist(g, names = FALSE),
          edge_attributes = igraph::edge_attr(g),
          graph_attributes = igraph::graph_attr(g)
        )
      })
    ))
  }
  if (is.data.frame(x)) {
    return(lapply(x, annotation_semantics))
  }
  x
}

expect_annotation_equal <- function(actual, expected) {
  expect_equal(
    annotation_semantics(suppressWarnings(actual)),
    annotation_semantics(suppressWarnings(expected))
  )
}

test_that("cache persists and unchanged versions do not rebuild", {
  local_cache_fixture()
  expect_identical(glyanno_cache_info()$status, "missing")
  expect_equal(
    as.character(comp_to_struc(c("H1N1", "Glc(2)"), return_best = TRUE)),
    c("Gal(b1-3)GalNAc(a1-", "Glc(a1-4)Glc(b1-")
  )
  expect_identical(file.exists(.annotation_cache_path()), FALSE)
  suppressMessages(build_glyanno_cache())
  expect_identical(glyanno_cache_info()$status, "ready")
  before <- readBin(
    .annotation_cache_path(),
    "raw",
    n = file.info(.annotation_cache_path())$size
  )
  .reset_annotation_cache()
  local_mocked_bindings(.build_annotation_view = function(...) {
    stop("must reuse disk cache")
  })
  expect_no_error(build_glyanno_cache())
  expect_equal(
    as.character(comp_to_struc(c("H1N1", "Glc(2)"), return_best = TRUE)),
    c("Gal(b1-3)GalNAc(a1-", "Glc(a1-4)Glc(b1-")
  )
  expect_identical(
    readBin(.annotation_cache_path(), "raw", n = length(before)),
    before
  )
  clear_glyanno_cache()
  expect_identical(glyanno_cache_info()$status, "missing")
  expect_identical(ls(.annotation_cache), character())
})

test_that("version changes invalidate disk and session caches", {
  local_cache_fixture()
  suppressMessages(build_glyanno_cache())
  .annotation_view("composition")
  metadata <- .cache_metadata()
  metadata$glydb <- "99.0.0"
  local_mocked_bindings(.cache_metadata = function() metadata)
  expect_identical(glyanno_cache_info()$status, "outdated")
  expect_null(.cache_session()$views[["composition_concrete"]])
  suppressMessages(build_glyanno_cache())
  expect_identical(glyanno_cache_info()$status, "ready")
  expect_identical(glyanno_cache_info()$stored$glydb, "99.0.0")
})

test_that("unreadable caches fall back and failed rebuilds preserve files", {
  local_cache_fixture()
  dir.create(dirname(.annotation_cache_path()))
  saveRDS(1, .annotation_cache_path())
  expect_identical(glyanno_cache_info()$status, "unreadable")
  expect_length(.annotation_view("composition")$db, 2L)
  suppressMessages(build_glyanno_cache())
  before <- tools::md5sum(.annotation_cache_path())
  local_mocked_bindings(.build_annotation_view = function(...) {
    stop("interrupted build")
  })
  expect_snapshot(
    error = TRUE,
    suppressMessages(build_glyanno_cache(force = TRUE))
  )
  expect_identical(tools::md5sum(.annotation_cache_path()), before)
  expect_identical(dir.exists(paste0(.annotation_cache_path(), ".lock")), FALSE)
})

test_that("cache startup reminder is actionable and suppressible", {
  local_cache_fixture()
  local_mocked_bindings(.cache_is_interactive = function() TRUE)
  expect_message(.onAttach(NULL, "glyanno"), "build_glyanno_cache")
  expect_no_message(suppressPackageStartupMessages(.onAttach(NULL, "glyanno")))
  suppressMessages(build_glyanno_cache())
  expect_no_message(.onAttach(NULL, "glyanno"))
  metadata <- .cache_metadata()
  metadata$glydb <- "99.0.0"
  local_mocked_bindings(.cache_metadata = function() metadata)
  expect_message(.onAttach(NULL, "glyanno"), "outdated.*99.0.0")
  local_mocked_bindings(.cache_is_interactive = function() FALSE)
  expect_no_message(.onAttach(NULL, "glyanno"))
})

test_that("build and clear respect a concurrent writer", {
  local_cache_fixture()
  dir.create(dirname(.annotation_cache_path()))
  lock <- paste0(.annotation_cache_path(), ".lock")
  dir.create(lock)
  expect_snapshot(error = TRUE, clear_glyanno_cache())
  expect_identical(dir.exists(lock), TRUE)
})

test_that("forced and outdated rebuilds work without overwriting rename destinations", {
  local_cache_fixture()
  local_mocked_bindings(.rename_cache_file = function(from, to) {
    if (file.exists(to)) {
      return(FALSE)
    }
    file.rename(from, to)
  })
  suppressMessages(build_glyanno_cache())
  expect_no_error(suppressMessages(build_glyanno_cache(force = TRUE)))
  expect_identical(glyanno_cache_info()$status, "ready")
  metadata <- .cache_metadata()
  metadata$glydb <- "99.0.0"
  local_mocked_bindings(.cache_metadata = function() metadata)
  expect_identical(glyanno_cache_info()$status, "outdated")
  expect_no_error(suppressMessages(build_glyanno_cache()))
  expect_identical(glyanno_cache_info()$stored, metadata)
  expect_identical(
    list.files(dirname(.annotation_cache_path()), pattern = "backup"),
    character()
  )
})

test_that("a failed replacement restores the previous cache", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  writeBin(charToRaw("previous cache"), "cache.rds")
  writeBin(charToRaw("replacement cache"), "new.rds")
  local_mocked_bindings(.rename_cache_file = function(from, to) {
    if (identical(from, "new.rds") || file.exists(to)) {
      return(FALSE)
    }
    file.rename(from, to)
  })
  expect_snapshot(error = TRUE, .replace_cache_file("new.rds", "cache.rds"))
  expect_identical(
    readBin("cache.rds", "raw", n = 100),
    charToRaw("previous cache")
  )
  expect_identical(list.files(pattern = "backup"), character())
})

test_that("metadata inspection and attachment never deserialize the payload", {
  local_cache_fixture()
  suppressMessages(build_glyanno_cache())
  local_mocked_bindings(
    .cache_is_interactive = function() TRUE,
    .read_cache_file = function(...) stop("payload must remain lazy")
  )
  expect_identical(glyanno_cache_info()$status, "ready")
  expect_no_message(.onAttach(NULL, "glyanno"))
  metadata <- .cache_metadata()
  metadata$glydb <- "99.0.0"
  local_mocked_bindings(.cache_metadata = function() metadata)
  expect_identical(glyanno_cache_info()$status, "outdated")
  expect_message(.onAttach(NULL, "glyanno"), "outdated")
})

test_that("payload loading occurs only once on the first annotation", {
  local_cache_fixture()
  suppressMessages(build_glyanno_cache())
  read <- .read_cache_file
  reads <- 0L
  local_mocked_bindings(.read_cache_file = function(path) {
    reads <<- reads + 1L
    read(path)
  })
  glyanno_cache_info()
  expect_identical(reads, 0L)
  comp_to_struc("H1N1")
  comp_to_struc("H1N1")
  expect_identical(reads, 1L)
})

test_that("truncated and legacy files are rejected without reading their payloads", {
  local_cache_fixture()
  suppressMessages(build_glyanno_cache())
  path <- .annotation_cache_path()
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  writeBin(head(bytes, -1L), path)
  expect_identical(glyanno_cache_info()$status, "unreadable")
  expect_null(.read_cache_file(path))
  saveRDS(list(metadata = .cache_metadata(), views = list()), path)
  expect_identical(glyanno_cache_info()$status, "unreadable")
  expect_null(.read_cache_file(path))
})

test_that("all persisted views preserve live database and annotation semantics", {
  withr::local_options(
    glyanno.cache_dir = tempfile("glyanno-parity-"),
    lifecycle_verbosity = "quiet"
  )
  withr::defer(unlink(getOption("glyanno.cache_dir"), recursive = TRUE))
  local_mocked_bindings(.annotation_cache = new.env(parent = emptyenv()))
  suppressMessages(build_glyanno_cache())
  .reset_annotation_cache()
  local_mocked_bindings(.build_annotation_view = function(...) {
    stop("unexpected rebuild")
  })
  for (kind in c("composition", "structure")) {
    levels <- if (kind == "structure") c("intact", "topological") else "intact"
    for (level in levels) {
      for (mono in c("concrete", "generic")) {
        args <- list(mono_type = mono)
        if (kind == "structure") {
          args$structure_level <- level
        }
        getter <- if (kind == "structure") {
          glydb::glydb_structures
        } else {
          glydb::glydb_compositions
        }
        for (filters in list(
          list(),
          list(species = "Homo sapiens"),
          list(glycan_type = "O"),
          list(species = "Mus musculus", glycan_type = "N"),
          list(mono_range = list(Hex = c(0L, Inf), HexNAc = c(0L, 6L))),
          list(mono_range = list(Hex = c(100L, 101L)))
        )) {
          filters <- c(args, filters)
          expected <- do.call(getter, filters)
          selected <- .cached_db(kind, filters)
          actual <- selected$view$db[selected$ids]
          expect_identical(as.character(actual), as.character(expected))
          expect_identical(
            attr(actual, "confidence"),
            attr(expected, "confidence")
          )
        }
      }
    }
  }
  filters <- list(species = "Homo sapiens", glycan_type = "O-GalNAc")
  comp_db <- do.call(glydb::glydb_compositions, filters)
  struc_db <- do.call(glydb::glydb_structures, filters)
  queries <- c("H1N1", "Gal(1)HexNAc(1)", "Gal(1)GalNAc(1)", NA, "H1N1")
  for (best in c(FALSE, TRUE)) {
    expect_annotation_equal(
      do.call(comp_to_struc, c(list(queries, return_best = best), filters)),
      comp_to_struc(queries, db = struc_db, return_best = best)
    )
    expect_equal(
      do.call(enhance_comp, c(list(queries, return_best = best), filters)),
      enhance_comp(queries, db = comp_db, return_best = best)
    )
    strucs <- c(
      "Hex(??-?)HexNAc(??-",
      NA,
      "Gal(??-?)GalNAc(??-",
      "Hex(??-?)HexNAc(??-"
    )
    expect_annotation_equal(
      do.call(enhance_struc, c(list(strucs, return_best = best), filters)),
      enhance_struc(strucs, db = struc_db, return_best = best)
    )
  }
  selected <- .cached_db("composition", filters)
  for (deriv in c("none", "permethyl", "peracetyl")) {
    for (type in c("mono", "average")) {
      dict <- glyanno_mass_dict(deriv, type)
      for (charge in c(-2L, 0L, 1L, 3L)) {
        adduct <- if (charge < 0L) "Cl-" else "Na+"
        expected <- suppressWarnings(calculate_mz(
          comp_db,
          charge,
          adduct,
          dict,
          safe = FALSE
        ))
        expect_equal(.cached_mz(selected, charge, adduct, dict), expected)
        expect_equal(
          do.call(
            mz_to_comp,
            c(
              list(
                expected[1:3],
                mass_dict = dict,
                charge = charge,
                adduct = adduct
              ),
              filters
            )
          ),
          mz_to_comp(
            expected[1:3],
            db = comp_db,
            mass_dict = dict,
            charge = charge,
            adduct = adduct
          )
        )
      }
    }
  }
  dict["Hex"] <- dict["Hex"] + 0.123
  expect_equal(
    .cached_mz(selected, 1L, "H+", dict),
    suppressWarnings(calculate_mz(comp_db, mass_dict = dict, safe = FALSE))
  )
  strucs <- glydb::glydb_data$glycan_structure[c(1L, 10L, 1L, NA_integer_)]
  expect_identical(
    local_struc_glytoucan_accessions(strucs),
    glydb::glydb_data$glytoucan_ac[match(
      strucs,
      glydb::glydb_data$glycan_structure
    )]
  )
})
