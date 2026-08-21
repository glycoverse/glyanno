#' Ensure the glycans is a [glyrepr::glycan_composition()]
#'
#' Convert all supported inputs into a [glyrepr::glycan_composition()].
#' @param glycans Glycans to process.
#' @param allow_structure Whether to allow structure inputs. Default is TRUE.
#' @returns A [glyrepr::glycan_composition()].
#' @noRd
.ensure_glycan_composition <- function(glycans, allow_structure = TRUE) {
  if (is.character(glycans)) {
    # ===== Character =====
    # Two cases: 1. Byonic composition strings, 2. structure strings
    # First assume composition strings
    tryCatch(
      return(glyrepr::as_glycan_composition(glycans)),
      error = function(e) NULL
    )
    # Then try to parse structure strings
    if (allow_structure) {
      tryCatch(
        {
          struc <- glyparse::auto_parse(glycans)
          return(glyrepr::as_glycan_composition(struc))
        },
        error = function(e) NULL
      )
    }
    # If both fail, raise an error
    if (allow_structure) {
      cli::cli_abort(
        "Cannot parse {.arg glycans} as glycan composition or structure strings."
      )
    } else {
      cli::cli_abort(
        "Cannot parse {.arg glycans} as glycan composition strings."
      )
    }
  } else if (glyrepr::is_glycan_composition(glycans)) {
    # ===== Glycan Composition =====
    glycans
  } else if (glyrepr::is_glycan_structure(glycans)) {
    # ===== Glycan Structure =====
    if (allow_structure) {
      glyrepr::as_glycan_composition(glycans)
    } else {
      cli::cli_abort("Cannot use glycan structures as input.")
    }
  } else {
    # ===== Other Types =====
    if (allow_structure) {
      cli::cli_abort(c(
        "{.arg glycans} must be a character vector, a {.fn glyrepr::glycan_composition} vector, or a {.fn glyrepr::glycan_structure} vector.",
        "x" = "Got {.cls {class(glycans)}}."
      ))
    } else {
      cli::cli_abort(c(
        "{.arg glycans} must be a character vector, a {.fn glyrepr::glycan_composition} vector.",
        "x" = "Got {.cls {class(glycans)}}."
      ))
    }
  }
}

#' Ensure the glycan structure is a [glyrepr::glycan_structure()]
#'
#' Convert all supported inputs into a [glyrepr::glycan_structure()].
#' @param strucs Glycan structures to process.
#' @returns A [glyrepr::glycan_structure()].
#' @noRd
.ensure_glycan_structure <- function(strucs) {
  if (is.character(strucs)) {
    tryCatch(
      return(glyparse::auto_parse(strucs)),
      error = function(e) {
        cli::cli_abort(
          "Cannot parse {.arg strucs} as glycan structure strings."
        )
      }
    )
  } else if (glyrepr::is_glycan_structure(strucs)) {
    strucs
  } else {
    cli::cli_abort(c(
      "{.arg strucs} must be a character vector or a {.fn glyrepr::glycan_structure} vector.",
      "x" = "Got {.cls {class(strucs)}}."
    ))
  }
}

#' Check if db has confidence attribute when return_best is TRUE
#' @noRd
.check_return_best_arg <- function(db, return_best, arg = "db") {
  if (isTRUE(return_best) && !.is_glydb_vector(db)) {
    cli::cli_abort(c(
      "{.arg {arg}} must have a {.val confidence} attribute when selecting the best match.",
      "i" = "Use {.fun glydb::glydb_compositions} or {.fun glydb::glydb_structures} to get a database with confidence scores."
    ))
  }
}

.prepare_struc_db <- function(db, arg = "db") {
  if (is.null(db)) {
    db <- glydb::glydb_structures(structure_level = "intact")
  } else {
    if (!.is_glydb_vector(db)) {
      if (length(db) == 0) {
        cli::cli_abort("{.arg db} cannot be of 0 length.")
      }
      db <- .ensure_glycan_structure(db)
      db <- unique(db)
    }
  }
  .drop_floating_structures(db, arg)
}

.prepare_denovo_struc_db <- function(fallback_db) {
  if (is.null(fallback_db)) {
    fallback_db <- glydb::glydb_structures(structure_level = "topological")
  }

  db <- .prepare_struc_db(fallback_db, arg = "fallback_db")
  mono_types <- glyrepr::get_mono_type(db)
  if (any(!is.na(mono_types) & mono_types != "concrete")) {
    cli::cli_abort(
      "{.arg fallback_db} must contain only concrete structures."
    )
  }

  confidence <- attr(db, "confidence")
  db <- glyrepr::remove_linkages(db)

  keys <- as.character(db)
  keep <- !duplicated(keys)
  if (!is.null(confidence)) {
    groups <- match(keys, keys[keep])
    confidence_groups <- split(
      confidence,
      factor(groups, levels = seq_along(keys[keep]))
    )
    confidence <- vapply(
      confidence_groups,
      function(values) {
        if (all(is.na(values))) {
          return(NA_real_)
        }
        max(values, na.rm = TRUE)
      },
      numeric(1)
    )
  }
  db <- db[keep]
  if (!is.null(confidence)) {
    attr(db, "confidence") <- unname(confidence)
  }
  db
}

.assert_denovo_strucs <- function(strucs) {
  present <- !is.na(strucs)
  mono_types <- glyrepr::get_mono_type(strucs)
  structure_levels <- glyrepr::get_structure_level(strucs)
  invalid <- present &
    (mono_types != "generic" | structure_levels != "topological")
  if (any(invalid)) {
    cli::cli_abort(c(
      "{.arg strucs} must contain only generic topological structures.",
      "i" = "Missing values are allowed, but every other element must have generic monosaccharides and no linkage information."
    ))
  }
  invisible(strucs)
}

.prepare_comp_db <- function(db) {
  if (is.null(db)) {
    db <- glydb::glydb_compositions()
  } else {
    if (!.is_glydb_vector(db)) {
      if (length(db) == 0) {
        cli::cli_abort("{.arg db} cannot be of 0 length.")
      }
      db <- .ensure_glycan_composition(db)
      db <- unique(db)
    }
  }
  db
}

.is_glydb_vector <- function(x) {
  !is.null(attr(x, "confidence"))
}

.prepare_result <- function(res, return_best, raw_col, new_col) {
  if (return_best) {
    res |>
      dplyr::arrange(.data$row_id, dplyr::desc(.data$confidence)) |>
      dplyr::group_by(.data$row_id) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::arrange(.data$row_id) |>
      dplyr::pull(.data[[new_col]])
  } else {
    res |>
      dplyr::filter(!is.na(.data[[new_col]])) |>
      dplyr::arrange(.data$row_id) |>
      dplyr::select(all_of(c(raw_col, new_col, "confidence")))
  }
}

.has_unresolved_floating <- function(strucs) {
  floating <- glyrepr::has_floating_parts(strucs) |
    glyrepr::has_floating_substituents(strucs)
  replace(floating, is.na(floating), FALSE)
}

.warn_floating_structures <- function(strucs, arg, action) {
  floating <- .has_unresolved_floating(strucs)
  if (any(floating)) {
    cli::cli_warn(c(
      "{.arg {arg}} contains {sum(floating)} structure{?s} with unresolved floating parts or substituents.",
      "i" = action
    ))
  }
  floating
}

.drop_floating_structures <- function(strucs, arg) {
  confidence <- attr(strucs, "confidence")
  floating <- .warn_floating_structures(
    strucs,
    arg,
    "Those database structures were excluded from matching."
  )
  strucs <- strucs[!floating]
  if (!is.null(confidence)) {
    attr(strucs, "confidence") <- confidence[!floating]
  }
  strucs
}

.replace_floating_structures <- function(strucs, arg = "strucs") {
  floating <- .warn_floating_structures(
    strucs,
    arg,
    paste0(
      "Those input structures were excluded from matching and are returned ",
      "as missing values in aligned outputs."
    )
  )
  if (!any(floating)) {
    return(list(strucs = strucs, rejected = floating))
  }
  iupacs <- as.character(strucs)
  iupacs[floating] <- NA_character_
  list(
    strucs = glyrepr::as_glycan_structure(iupacs),
    rejected = floating
  )
}

.new_composition_match_index <- function(candidates) {
  list(
    generic_keys = as.character(glyrepr::convert_to_generic(candidates)),
    counts = as.list(candidates)
  )
}

.composition_match_ids <- function(pattern, index) {
  if (length(pattern) == 0 || is.na(pattern)) {
    return(integer())
  }

  generic_key <- as.character(glyrepr::convert_to_generic(pattern))
  candidate_ids <- which(index$generic_keys == generic_key)
  if (length(candidate_ids) == 0) {
    return(integer())
  }

  pattern_counts <- as.list(pattern)[[1]]
  mono_types <- glyrepr::get_mono_type(names(pattern_counts))
  concrete_counts <- pattern_counts[mono_types == "concrete"]
  if (length(concrete_counts) == 0) {
    return(candidate_ids)
  }

  compatible <- vapply(
    index$counts[candidate_ids],
    function(candidate_counts) {
      all(vapply(
        names(concrete_counts),
        function(mono) {
          candidate_count <- unname(candidate_counts[mono])
          if (is.na(candidate_count)) {
            candidate_count <- 0L
          }
          candidate_count >= concrete_counts[[mono]]
        },
        logical(1)
      ))
    },
    logical(1)
  )
  candidate_ids[compatible]
}

.best_match_id <- function(candidate_ids, best_rank) {
  if (length(candidate_ids) == 0) {
    return(NA_integer_)
  }
  candidate_ids[[which.min(best_rank[candidate_ids])]]
}

.assert_concrete <- function(strucs) {
  if (length(strucs) == 0) {
    return(NULL)
  }
  mono_types <- glyrepr::get_mono_type(strucs)
  if (any(!is.na(mono_types) & mono_types != "concrete")) {
    cli::cli_abort(
      "{.arg strucs} must have concrete monosaccharides (e.g. Gal, GalNAc)."
    )
  }
}
