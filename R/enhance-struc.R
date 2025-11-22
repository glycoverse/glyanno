#' Enhance glycan structure
#'
#' Given a glycan structure of any resolution level (see [glyrepr::get_structure_level()] for details),
#' this function gives all possible glycan structures of higher resolution level.
#'
#' @param strucs A [glyrepr::glycan_structure()] vector,
#'   or a character vector of glycan structure strings supported by [glyparse::auto_parse()].
#'   Glycan structures with level higher or same as `to_level` will be returned as is.
#'   Glycan structures with level lower than `to_level` will be enhanced to the level of `to_level`.
#' @param to_level The resolution level to enhance to.
#'   Can be "intact" or "topological". Default is "intact".
#' @param db A [glydb::glydb_structures()] vector,
#'   or a character vector of glycan structure strings supported by [glyparse::auto_parse()].
#'   All structures in `db` must be at the same resolution level as `to_level`.
#'   If not provided, a default structure vector is loaded from [glydb::glydb_structures()]
#'   with the specified `to_level`.
#'
#' @returns A tibble with the following columns:
#'   - `raw`: The original glycan structures.
#'   - `enhanced`: The enhanced glycan structures.
#'     Note that one `raw` glycan structure can have different `enhanced` glycan structures
#'     as multiple rows in the result.
#'
#' @examples
#' # From topological level to intact level
#' enhance_struc("Gal(??-?)GalNAc(??-", to_level = "intact")
#'
#' # From basic level to topological level
#' enhance_struc("Hex(??-?)HexNAc(??-", to_level = "topological")
#'
#' # From partial level to intact level
#' enhance_struc("Gal(b1-?)GalNAc(a1-", to_level = "intact")
#'
#' @export
enhance_struc <- function(strucs, to_level = "intact", db = NULL) {
  # Input validation and preparation
  strucs <- .ensure_glycan_structure(strucs)
  checkmate::assert_choice(to_level, c("intact", "topological"))
  if (is.null(db)) {
    db <- glydb::glydb_structures(structure_level = to_level)
  } else {
    db <- .ensure_glycan_structure(db)
  }
  db <- unique(db)
  db_struc_level <- glyrepr::get_structure_level(db)
  if (any(db_struc_level != to_level)) {
    cli::cli_warn("Some structures in `db` are not at the same resolution level as `to_level`, which will be dropped.")
    db <- db[db_struc_level == to_level]
  }

  # Define level ranks (higher number means higher resolution)
  level_ranks <- c("basic" = 1, "topological" = 2, "partial" = 3, "intact" = 4)

  target_rank <- level_ranks[[to_level]]
  struc_levels <- glyrepr::get_structure_level(strucs)
  struc_ranks <- level_ranks[struc_levels]

  # Create a dataframe to track original IDs
  strucs_df <- tibble::tibble(
    raw = strucs,
    row_id = seq_along(strucs),
    rank = struc_ranks
  )

  # Identify structures to enhance
  keep_df <- strucs_df |> dplyr::filter(.data$rank >= target_rank)
  enhance_df <- strucs_df |> dplyr::filter(.data$rank < target_rank)

  res_keep <- NULL
  res_enhanced <- NULL

  # Process structures to keep
  if (nrow(keep_df) > 0) {
    res_keep <- keep_df |>
      dplyr::select(all_of(c("raw", "row_id"))) |>
      dplyr::mutate(enhanced = .data$raw)
  }

  # Process structures to enhance
  if (nrow(enhance_df) > 0) {
    to_enhance <- enhance_df$raw

    if (length(db) == 0) {
      # Empty db, no matches
      res_enhanced <- tibble::tibble(
        raw = to_enhance[integer(0)],
        enhanced = db[integer(0)],
        row_id = integer(0)
      )
    } else {
      # Check matches
      # db are targets (glycans), to_enhance are patterns (motifs)
      matches <- glymotif::have_motifs(db, to_enhance, alignments = "whole")

      # matches: rows=db, cols=to_enhance
      match_indices <- which(matches, arr.ind = TRUE)

      if (nrow(match_indices) > 0) {
        # match_indices[, "col"] are indices into to_enhance
        # match_indices[, "row"] are indices into db

        # Get row_ids from enhance_df
        matched_row_ids <- enhance_df$row_id[match_indices[, "col"]]
        matched_raw <- enhance_df$raw[match_indices[, "col"]]
        matched_enhanced <- db[match_indices[, "row"]]

        res_enhanced <- tibble::tibble(
          raw = matched_raw,
          enhanced = matched_enhanced,
          row_id = matched_row_ids
        )
      } else {
        # No matches found
        res_enhanced <- tibble::tibble(
          raw = to_enhance[integer(0)],
          enhanced = db[integer(0)],
          row_id = integer(0)
        )
      }
    }
  }

  # Combine results
  res <- dplyr::bind_rows(res_keep, res_enhanced) |>
    dplyr::arrange(.data$row_id) |>
    dplyr::select(all_of(c("raw", "enhanced")))

  # Handle empty result
  if (nrow(res) == 0) {
    res <- tibble::tibble(
      raw = glyrepr::glycan_structure(),
      enhanced = glyrepr::glycan_structure()
    )
  }

  res
}
