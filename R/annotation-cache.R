.view_key <- function(kind, mono_type, structure_level = "intact") {
  paste(
    c(kind, if (kind == "structure") structure_level, mono_type),
    collapse = "_"
  )
}

.annotation_view <- function(
  kind,
  mono_type = "concrete",
  structure_level = "intact"
) {
  state <- .cache_session()
  key <- .view_key(kind, mono_type, structure_level)
  if (is.null(state$views[[key]])) {
    state$views[[key]] <- .build_annotation_view(
      kind,
      mono_type,
      structure_level
    )
  }
  state$views[[key]]
}

.build_annotation_view <- function(kind, mono_type, structure_level) {
  getter <- if (kind == "structure") {
    glydb::glydb_structures
  } else {
    glydb::glydb_compositions
  }
  args <- list(mono_type = mono_type)
  if (kind == "structure") {
    args$structure_level <- structure_level
  }
  db <- do.call(getter, args)
  keys <- as.character(db)
  composition <- if (kind == "structure") {
    glyrepr::as_glycan_composition(db)
  } else {
    db
  }
  generic <- glyrepr::convert_to_generic(composition)
  membership <- function(field, choices) {
    stats::setNames(
      lapply(choices, function(choice) {
        selected <- do.call(
          getter,
          c(args, stats::setNames(list(choice), field))
        )
        match(as.character(selected), keys)
      }),
      choices
    )
  }
  confidence <- attr(db, "confidence")
  counts <- as.list(composition)
  generic_counts <- as.list(generic)
  mass_counts <- .mass_count_matrix(generic)
  dictionaries <- expand.grid(
    deriv = c("none", "permethyl", "peracetyl"),
    type = c("mono", "average"),
    stringsAsFactors = FALSE
  )
  masses <- lapply(seq_len(nrow(dictionaries)), function(i) {
    dict <- glyanno_mass_dict(dictionaries$deriv[i], dictionaries$type[i])
    list(dictionary = dict, residue = .residue_mass(mass_counts, dict))
  })
  list(
    mono_types = glyrepr::get_mono_type(composition),
    db = db,
    composition = composition,
    exact = .new_exact_composition_match_index(composition),
    compatible = .new_composition_match_index(
      composition,
      as.character(generic)
    ),
    best_rank = match(
      seq_along(db),
      order(
        replace(confidence, is.na(confidence), -Inf),
        decreasing = TRUE,
        method = "radix"
      )
    ),
    floating = if (kind == "structure") {
      .has_unresolved_floating(db)
    } else {
      rep(FALSE, length(db))
    },
    counts = counts,
    generic_counts = generic_counts,
    mass_counts = mass_counts,
    species = membership("species", glydb::glydb_species()),
    types = membership(
      "glycan_type",
      c(
        "HMO",
        "N",
        "GSL",
        "GAG",
        "O",
        "O-GalNAc",
        "O-GlcNAc",
        "O-Man",
        "O-Fuc",
        "O-Glc",
        "GPI"
      )
    ),
    masses = masses
  )
}

.mass_count_matrix <- function(composition) {
  components <- c(
    glyrepr::available_monosaccharides("generic"),
    glyrepr::available_substituents()
  )
  result <- vapply(
    components,
    function(component) {
      as.numeric(glyrepr::count_mono(composition, component))
    },
    numeric(length(composition))
  )
  matrix(
    result,
    nrow = length(composition),
    ncol = length(components),
    dimnames = list(NULL, components)
  )
}

.residue_mass <- function(counts, dictionary) {
  present <- colnames(counts) %in% names(dictionary)
  # Keep calculate_mz's summation order, including zero-count components.
  values <- sweep(
    counts[, present, drop = FALSE],
    2,
    dictionary[colnames(counts)[present]],
    `*`
  )
  mass <- rowSums(values)
  mass[rowSums(counts[, !present, drop = FALSE], na.rm = TRUE) > 0] <- NA_real_
  unname(mass)
}

.build_accession_lookup <- function() {
  keys <- as.character(glydb::glydb_data$glycan_structure)
  keep <- !duplicated(keys)
  stats::setNames(glydb::glydb_data$glytoucan_ac[keep], keys[keep])
}

.cached_accessions <- function(strucs) {
  state <- .cache_session()
  if (is.null(state$accessions)) {
    state$accessions <- .build_accession_lookup()
  }
  unname(state$accessions[match(as.character(strucs), names(state$accessions))])
}

.cached_db <- function(kind, filters = list()) {
  mono <- if ("mono_type" %in% names(filters)) filters$mono_type else "concrete"
  level <- if ("structure_level" %in% names(filters)) {
    filters$structure_level
  } else {
    "intact"
  }
  checkmate::assert_choice(mono, c("concrete", "generic"))
  checkmate::assert_choice(level, c("intact", "topological"))
  view <- .annotation_view(kind, mono, level)
  checkmate::assert_choice(filters$species, names(view$species), null.ok = TRUE)
  checkmate::assert_choice(
    filters$glycan_type,
    names(view$types),
    null.ok = TRUE
  )
  ids <- seq_along(view$db)
  if (!is.null(filters$species)) {
    ids <- ids[ids %in% view$species[[filters$species]]]
  }
  if (!is.null(filters$glycan_type)) {
    ids <- ids[ids %in% view$types[[filters$glycan_type]]]
  }
  if (!is.null(filters$mono_range)) {
    .validate_cached_mono_range(filters$mono_range, mono)
    counts <- if (
      glyrepr::get_mono_type(names(filters$mono_range)[1]) == "generic"
    ) {
      view$generic_counts
    } else {
      view$counts
    }
    ranges <- filters$mono_range
    keep <- vapply(
      counts[ids],
      function(x) {
        if (length(setdiff(names(x), names(ranges)))) {
          return(FALSE)
        }
        all(vapply(
          names(ranges),
          function(mono) {
            n <- unname(x[mono])
            if (is.na(n)) {
              n <- 0L
            }
            n >= ranges[[mono]][1] && n <= ranges[[mono]][2]
          },
          logical(1)
        ))
      },
      logical(1)
    )
    ids <- ids[keep]
  }
  structure(list(view = view, ids = ids), class = "glyanno_cached_db")
}

.is_cached_db <- function(db) inherits(db, "glyanno_cached_db")

.cached_structure_ids <- function(db) {
  db$ids[!db$view$floating[db$ids]]
}

.cached_match_index <- function(db, mode = "compatible", kind = "composition") {
  ids <- if (kind == "structure") .cached_structure_ids(db) else db$ids
  index <- db$view[[mode]]
  # Store full-database IDs; map only the query's matches into the selected view.
  index$active_ids <- ids
  list(
    db = .subset_cached_db(db$view$db, ids),
    match_index = index,
    best_rank = db$view$best_rank[ids]
  )
}

.cached_mz <- function(db, charge, adduct, mass_dict) {
  view <- db$view
  found <- which(vapply(
    view$masses,
    function(x) identical(x$dictionary, mass_dict),
    logical(1)
  ))
  residue <- if (length(found)) {
    view$masses[[found[1]]]$residue
  } else {
    .residue_mass(view$mass_counts, mass_dict)
  }
  mz <- residue[db$ids] + mass_dict[adduct] * abs(charge) + mass_dict["red_end"]
  if (charge != 0) {
    mz <- mz / abs(charge)
  }
  unname(mz)
}

# Mirrors the public glydb mono_range contract.
.validate_cached_mono_range <- function(mono_range, mono_type) {
  # NULL is valid (no filtering)
  if (is.null(mono_range)) {
    return(invisible(NULL))
  }

  # Check that mono_range is a non-empty list
  checkmate::assert_list(mono_range)
  if (length(mono_range) == 0L) {
    cli::cli_abort(c(
      "{.arg mono_range} must contain at least one named monosaccharide range.",
      "i" = "Provide one or more named entries, e.g. {.code list(Hex = c(3L, 9L))}, or use {.code NULL} for no filtering."
    ))
  }

  # Check for duplicate names
  names_list <- names(mono_range)
  if (is.null(names_list) || any(names_list == "")) {
    cli::cli_abort(c(
      "All elements in {.arg mono_range} must be named.",
      "i" = "Use named list entries like {.code list(Hex = c(3L, 9L))}."
    ))
  }

  dup_names <- names_list[duplicated(names_list)]
  if (length(dup_names) > 0) {
    cli::cli_abort(c(
      "Monosaccharide names in {.arg mono_range} must not be duplicated.",
      "x" = "Duplicated name{?s}: {.val {unique(dup_names)}}."
    ))
  }

  # Check that all names are valid monosaccharides
  valid_monos <- glyrepr::available_monosaccharides()
  invalid_monos <- setdiff(names_list, valid_monos)
  if (length(invalid_monos) > 0) {
    cli::cli_abort(c(
      "All monosaccharide names in {.arg mono_range} must be known.",
      "x" = "Unknown name{?s}: {.val {invalid_monos}}."
    ))
  }

  # Check that all names are of the same mono_type
  range_mono_types <- glyrepr::get_mono_type(names_list)
  if (length(unique(range_mono_types)) != 1) {
    cli::cli_abort(c(
      "All monosaccharide names in {.arg mono_range} must be of the same type (generic or concrete).",
      "x" = "Found both types."
    ))
  }

  # Special case: When `mono_type` is generic, mono type of the range names must be generic.
  range_mono_type <- glyrepr::get_mono_type(names(mono_range)[[1]])
  if (mono_type == "generic" && range_mono_type == "concrete") {
    cli::cli_abort(c(
      "Monosaccharide names in {.arg mono_range} must be {.val generic} (e.g. 'Hex') when {.arg mono_type} is {.val generic}."
    ))
  }

  # Validate each range element
  for (mono_name in names_list) {
    range_val <- mono_range[[mono_name]]

    # Check that element is a numeric vector of length 2
    if (!checkmate::test_numeric(range_val, len = 2)) {
      cli::cli_abort(c(
        "Each element in {.arg mono_range} must be a numeric vector of length 2.",
        "x" = "{.val {mono_name}} is not a valid numeric vector of length 2."
      ))
    }

    min_val <- range_val[1]
    max_val <- range_val[2]

    # Check that min is an integer (not double, not Inf)
    # Special case: c(0L, Inf) creates a double vector due to R's type coercion
    # We allow this when max is Inf and min is a whole number
    is_special_case <- is.infinite(max_val) &&
      is.double(min_val) &&
      min_val == as.integer(min_val)

    if (is.infinite(min_val)) {
      cli::cli_abort(c(
        "Minimum value in {.arg mono_range} must be an integer.",
        "x" = "{.val {mono_name}} has non-integer minimum value {.val {min_val}}.",
        "i" = "Use integer literals like {.code 3L} instead of {.code 3}."
      ))
    }

    if (!is.integer(min_val) && !is_special_case) {
      cli::cli_abort(c(
        "Minimum value in {.arg mono_range} must be an integer.",
        "x" = "{.val {mono_name}} has non-integer minimum value {.val {min_val}}.",
        "i" = "Use integer literals like {.code 3L} instead of {.code 3}."
      ))
    }

    # Check that max is either an integer or Inf
    is_valid_max <- is.infinite(max_val) || is.integer(max_val)
    if (!is_valid_max) {
      cli::cli_abort(c(
        "Maximum value in {.arg mono_range} must be an integer or {.val {Inf}}.",
        "x" = "{.val {mono_name}} has invalid maximum value {.val {max_val}}.",
        "i" = "Use integer literals like {.code 3L} instead of {.code 3}."
      ))
    }

    # Check min constraints: must be >= 0
    if (min_val < 0L) {
      cli::cli_abort(c(
        "Minimum value in {.arg mono_range} must be non-negative.",
        "x" = "{.val {mono_name}} has minimum value {.val {min_val}}."
      ))
    }

    # Check max constraints: must be >= 0 (Inf is allowed)
    if (max_val < 0L) {
      cli::cli_abort(c(
        "Maximum value in {.arg mono_range} must be non-negative.",
        "x" = "{.val {mono_name}} has maximum value {.val {max_val}}."
      ))
    }

    # Check that max >= min
    if (max_val < min_val) {
      cli::cli_abort(c(
        "Maximum value must be greater than or equal to minimum value in {.arg mono_range}.",
        "x" = "{.val {mono_name}} has max ({.val {max_val}}) < min ({.val {min_val}})."
      ))
    }
  }

  invisible(NULL)
}

.subset_cached_db <- function(db, ids) {
  result <- db[ids]
  attr(result, "confidence") <- attr(db, "confidence")[ids]
  result
}
