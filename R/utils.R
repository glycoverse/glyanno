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
.check_confidence_attr <- function(db, return_best) {
  if (isTRUE(return_best) && is.null(attr(db, "confidence"))) {
    cli::cli_abort(c(
      "`db` must have a {.val confidence} attribute when {.val return_best} is {.val TRUE}.",
      "i" = "Use {.fun glydb::glydb_compositions} or {.fun glydb::glydb_structures} to get a database with confidence scores."
    ))
  }
}
