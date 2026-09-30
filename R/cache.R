#' Manage the persistent annotation cache
#'
#' @description
#' `build_glyanno_cache()` is the magic that makes `glyanno` fast.
#' It precalculates database views and annotation indexes from the installed `glydb`.
#' You will only need to call this function once on your computer.
#' When `glydb` is updated, you'll also need to call this function to update the cache.
#' But don't worry. We'll remind you then.
#'
#' - `build_glyanno_cache()` reuses a valid cache unless `force = TRUE`.
#' - `glyanno_cache_info()` inspects it without building it, and
#' - `clear_glyanno_cache()` removes it and clears the session cache.
#'
#' @param force Rebuild even when the cache is current.
#' @details
#' The cache lives in `tools::R_user_dir("glyanno", "cache")`. Set the
#' `glyanno.cache_dir` option to use another directory. Only an explicit call
#' to `build_glyanno_cache()` writes persistent data. Without a disk cache,
#' annotation builds the required views in memory for the current session.
#'
#' On interactive attachment, a startup message reminds you to build a missing,
#' outdated, or unreadable cache. No computation or questions block attachment.
#' The cache is invalidated by a changed glydb version, glyrepr or igraph version,
#' R major/minor version, collation locale, or internal cache schema. Restart R
#' after updating dependencies. Development data changes without a version bump
#' require `force = TRUE`.
#'
#' Cached data include concrete/generic compositions, intact/topological
#' structures, confidence ranks, exact and compatible composition indexes,
#' residue counts, species/type membership, floating-structure flags, masses for
#' all built-in dictionaries, and local GlyTouCan accession keys. Structure
#' vectors retain their parsed graphs. Arbitrary structure matching still uses
#' glymotif at query time. Custom databases are never persisted.
#'
#' A complete replacement is written before the previous cache is replaced;
#' failed builds preserve the previous file. Only one completed cache is kept.
#'
#' @returns `build_glyanno_cache()` invisibly returns the cache file path.
#'   `glyanno_cache_info()` returns a list with `path`, `status` (one of
#'   `"missing"`, `"outdated"`, `"unreadable"`, or `"ready"`), `size_bytes`,
#'   `expected`, and `stored` version metadata.
#'   `clear_glyanno_cache()` invisibly returns `NULL`.
#' @examples
#' glyanno_cache_info()
#' @export
build_glyanno_cache <- function(force = FALSE) {
  checkmate::assert_flag(force)
  info <- glyanno_cache_info()
  if (!force && info$status == "ready") {
    return(invisible(info$path))
  }
  path <- info$path
  dir <- dirname(path)
  if (!dir.exists(dir) && !dir.create(dir, recursive = TRUE)) {
    cli::cli_abort(
      "Cannot create the {.pkg glyanno} cache directory {.file {dir}}."
    )
  }
  # A directory lock prevents concurrent writers from replacing one another.
  lock <- paste0(path, ".lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    cli::cli_abort(c(
      "Another cache build may be running.",
      "i" = "If it was interrupted, remove {.file {lock}} before retrying."
    ))
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  payload <- list(metadata = .cache_metadata(), views = list())
  for (kind in c("composition", "structure")) {
    levels <- if (kind == "structure") c("intact", "topological") else "intact"
    for (level in levels) {
      for (mono in c("concrete", "generic")) {
        key <- .view_key(kind, mono, level)
        cli::cli_inform(c("i" = "Preparing {.val {key}}..."))
        payload$views[[key]] <- .build_annotation_view(kind, mono, level)
      }
    }
  }
  payload$accessions <- .build_accession_lookup()
  tmp <- tempfile("annotation-", tmpdir = dir)
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(payload, tmp, compress = TRUE)
  if (!.valid_cache_payload(readRDS(tmp))) {
    cli::cli_abort("Cache validation failed; the previous cache is unchanged.")
  }
  if (!file.rename(tmp, path)) {
    cli::cli_abort("Cannot replace the cache file {.file {path}}.")
  }
  .reset_annotation_cache()
  cli::cli_inform(c("v" = "Annotation cache saved to {.file {path}}."))
  invisible(path)
}

#' @rdname build_glyanno_cache
#' @export
glyanno_cache_info <- function() {
  path <- .annotation_cache_path()
  expected <- .cache_metadata()
  stored <- NULL
  status <- "missing"
  if (file.exists(path)) {
    payload <- .read_cache_file(path)
    stored <- if (is.list(payload)) payload$metadata else NULL
    status <- if (!.valid_cache_payload(payload)) {
      "unreadable"
    } else if (!identical(stored, expected)) {
      "outdated"
    } else {
      "ready"
    }
  }
  list(
    path = path,
    status = status,
    size_bytes = if (file.exists(path)) unname(file.info(path)$size) else 0,
    expected = expected,
    stored = stored
  )
}

#' @rdname build_glyanno_cache
#' @export
clear_glyanno_cache <- function() {
  path <- .annotation_cache_path()
  if (dir.exists(paste0(path, ".lock"))) {
    cli::cli_abort("Cannot clear the cache while a build lock exists.")
  }
  if (file.exists(path) && unlink(path) != 0L) {
    cli::cli_abort("Cannot remove the cache file {.file {path}}.")
  }
  .reset_annotation_cache()
  invisible(NULL)
}

.annotation_cache <- new.env(parent = emptyenv())

.annotation_cache_path <- function() {
  dir <- getOption("glyanno.cache_dir", tools::R_user_dir("glyanno", "cache"))
  checkmate::assert_string(dir, min.chars = 1L)
  file.path(path.expand(dir), "annotation-v1.rds")
}

# Bump schema when cached representations or mass dictionary semantics change.
.cache_metadata <- function() {
  list(
    schema = 1L,
    glydb = as.character(getNamespaceVersion("glydb")),
    glyrepr = as.character(getNamespaceVersion("glyrepr")),
    igraph = as.character(getNamespaceVersion("igraph")),
    R = paste(
      R.version$major,
      strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1],
      sep = "."
    ),
    collation = Sys.getlocale("LC_COLLATE")
  )
}

.read_cache_file <- function(path) {
  tryCatch(suppressWarnings(readRDS(path)), error = function(e) NULL)
}

.valid_cache_payload <- function(x) {
  is.list(x) &&
    is.list(x$metadata) &&
    is.list(x$views) &&
    all(
      c(
        "composition_concrete",
        "composition_generic",
        "structure_intact_concrete",
        "structure_intact_generic",
        "structure_topological_concrete",
        "structure_topological_generic"
      ) %in%
        names(x$views)
    ) &&
    all(vapply(
      x$views,
      function(v) {
        is.list(v) &&
          all(
            c(
              "db",
              "composition",
              "mono_types",
              "mass_counts",
              "exact",
              "compatible",
              "best_rank",
              "floating",
              "counts",
              "generic_counts",
              "species",
              "types",
              "masses"
            ) %in%
              names(v)
          )
      },
      logical(1)
    )) &&
    is.character(x$accessions)
}

.reset_annotation_cache <- function() {
  rm(list = ls(.annotation_cache, all.names = TRUE), envir = .annotation_cache)
}

.cache_session <- function() {
  identity <- list(
    path = .annotation_cache_path(),
    metadata = .cache_metadata()
  )
  if (!identical(.annotation_cache$identity, identity)) {
    .reset_annotation_cache()
    .annotation_cache$identity <- identity
    payload <- .read_cache_file(identity$path)
    if (
      .valid_cache_payload(payload) &&
        identical(payload$metadata, identity$metadata)
    ) {
      .annotation_cache$views <- payload$views
      .annotation_cache$accessions <- payload$accessions
    } else {
      .annotation_cache$views <- list()
    }
  }
  .annotation_cache
}

.cache_is_interactive <- function() interactive()

.onAttach <- function(libname, pkgname) {
  if (!.cache_is_interactive()) {
    return(invisible(NULL))
  }
  tryCatch(
    {
      info <- glyanno_cache_info()
      if (info$status != "ready") {
        cli::cli_inform(
          c(
            "i" = "{.pkg glyanno}: annotation cache is {info$status} for {.pkg glydb} {info$expected$glydb}.",
            "i" = "Run {.code glyanno::build_glyanno_cache()} to build it."
          ),
          class = "packageStartupMessage"
        )
      }
    },
    error = function(e) {
      cli::cli_inform(
        c(
          "i" = "{.pkg glyanno}: cannot inspect the annotation cache.",
          "i" = "Run {.code glyanno::glyanno_cache_info()} for details."
        ),
        class = "packageStartupMessage"
      )
    }
  )
}
