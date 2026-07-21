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
.check_return_best_arg <- function(db, return_best) {
  if (isTRUE(return_best) && !.is_glydb_vector(db)) {
    cli::cli_abort(c(
      "`db` must have a {.val confidence} attribute when {.val return_best} is {.val TRUE}.",
      "i" = "Use {.fun glydb::glydb_compositions} or {.fun glydb::glydb_structures} to get a database with confidence scores."
    ))
  }
}

.prepare_struc_db <- function(db) {
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
  db
}

.prepare_denovo_struc_db <- function(db) {
  if (is.null(db)) {
    return(glydb::glydb_structures(structure_level = "topological"))
  }

  db <- .prepare_struc_db(db)
  db_level <- glyrepr::get_structure_level(db)
  if (identical(db_level, "basic")) {
    cli::cli_abort(
      "{.arg db} cannot contain basic structures when {.code method = \"denovo\"}."
    )
  }

  confidence <- attr(db, "confidence")
  if (db_level %in% c("partial", "intact")) {
    db <- glyrepr::reduce_structure_level(db, "topological")
  }

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
      dplyr::select(all_of(c(raw_col, new_col)))
  }
}

.assert_concrete <- function(strucs) {
  if (length(strucs) == 0) {
    return(NULL)
  }
  if (glyrepr::get_mono_type(strucs) != "concrete") {
    cli::cli_abort(
      "{.arg strucs} must have concrete monosaccharides (e.g. Gal, GalNAc)."
    )
  }
}
