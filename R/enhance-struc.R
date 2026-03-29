#' Enhance glycan structure
#'
#' Given a glycan structure of any resolution level (see [glyrepr::get_structure_level()] for details),
#' this function gives all possible glycan structures of higher resolution level.
#'
#' The target resolution level is determined from `db`. All structures in `db` must be at the same level.
#' When `db` is NULL, the default [glydb::glydb_structures()] at "intact" level is used.
#'
#' @param strucs A [glyrepr::glycan_structure()] vector,
#'   or a character vector of glycan structure strings supported by [glyparse::auto_parse()].
#'   Glycan structures with level higher or same as the level of `db` will be returned as is.
#'   Glycan structures with level lower than the level of `db` will be enhanced to that level.
#' @param db A [glydb::glydb_structures()] vector,
#'   or a character vector of glycan structure strings supported by [glyparse::auto_parse()].
#'   All structures in `db` must be at the same resolution level.
#'   If not provided, a default structure vector is loaded from [glydb::glydb_structures()]
#'   at "intact" level.
#' @param return_best Logical. If `TRUE`, only return the best matching structure
#'   (highest confidence) for each input structure. Requires `db` to have a
#'   `confidence` attribute. Use [glydb::glydb_structures()] for `db` to enable this feature.
#'   Default is `FALSE`.
#'
#' @returns If `return_best=TRUE`:
#'   An unnamed [glyrepr::glycan_structure()] vector with the same length as `strucs`.
#'   Unmatched structures are returned as `NA`.
#'   If `return_best=FALSE`:
#'   A tibble with the following columns:
#'   - `raw`: The original glycan structures.
#'   - `enhanced`: The enhanced glycan structures.
#'     Note that one `raw` glycan structure can have different `enhanced` glycan structures
#'     as multiple rows in the result.
#'
#' @examples
#' # From topological level to intact level
#' db_intact <- c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
#' enhance_struc("Gal(??-?)GalNAc(??-", db = db_intact)
#'
#' # From basic level to topological level
#' db_topo <- "Gal(??-?)GalNAc(??-"
#' enhance_struc("Hex(??-?)HexNAc(??-", db = db_topo)
#'
#' # From partial level to intact level
#' enhance_struc("Gal(b1-?)GalNAc(a1-", db = db_intact)
#'
#' @export
enhance_struc <- function(
  strucs,
  db = NULL,
  return_best = FALSE
) {
  # Input validation and preparation
  strucs <- .ensure_glycan_structure(strucs)
  checkmate::assert_flag(return_best)

  if (is.null(db)) {
    db <- glydb::glydb_structures(structure_level = "intact")
    confidences <- attr(db, "confidence")
  } else {
    confidences <- attr(db, "confidence")
    if (is.null(confidences)) {
      db <- .ensure_glycan_structure(db)
      db <- unique(db)
    }
  }

  .check_confidence_attr(confidences, return_best)

  # Determine target level from db
  db_struc_levels <- glyrepr::get_structure_level(db)
  unique_levels <- unique(db_struc_levels)
  if (length(unique_levels) > 1) {
    cli::cli_abort(c(
      "All structures in `db` must have the same structure level.",
      "x" = "Found {length(unique_levels)} different levels: {.val {unique_levels}}."
    ))
  }
  to_level <- unique_levels[1]

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
      # Empty db, all unmatched
      if (return_best) {
        na_vec <- rep(glyrepr::glycan_structure(NA), length(to_enhance))
        res_enhanced <- tibble::tibble(
          raw = to_enhance,
          enhanced = na_vec,
          row_id = enhance_df$row_id
        )
      } else {
        res_enhanced <- tibble::tibble(
          raw = to_enhance[integer(0)],
          enhanced = db[integer(0)],
          row_id = integer(0)
        )
      }
    } else {
      # Check matches
      # db are targets (glycans), to_enhance are patterns (motifs)
      matches <- glymotif::have_motifs(db, to_enhance, alignments = "whole")

      # For return_best=TRUE: one row per to_enhance, best match or NA
      # For return_best=FALSE: one row per match
      if (return_best) {
        # Build result for each pattern - one best match or NA
        res_list <- lapply(seq_along(to_enhance), function(i) {
          col_matches <- which(matches[, i])
          if (length(col_matches) > 0) {
            # Find best match by confidence
            confs <- confidences[col_matches]
            best_col <- col_matches[which.max(confs)]
            tibble::tibble(
              raw = to_enhance[i],
              enhanced = db[best_col],
              row_id = enhance_df$row_id[i]
            )
          } else {
            # No match - NA
            tibble::tibble(
              raw = to_enhance[i],
              enhanced = glyrepr::glycan_structure(NA),
              row_id = enhance_df$row_id[i]
            )
          }
        })
        res_enhanced <- dplyr::bind_rows(res_list)
      } else {
        # matches: rows=db, cols=to_enhance
        match_indices <- which(matches, arr.ind = TRUE)

        if (nrow(match_indices) > 0) {
          # match_indices[, "col"] are indices into to_enhance
          # match_indices[, "row"] are indices into db
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
  }

  # Combine results
  res <- dplyr::bind_rows(res_keep, res_enhanced) |>
    dplyr::arrange(.data$row_id) |>
    dplyr::select(all_of(c("raw", "enhanced", "row_id")))

  # Handle empty result
  if (nrow(res) == 0) {
    res <- tibble::tibble(
      raw = glyrepr::glycan_structure(),
      enhanced = glyrepr::glycan_structure(),
      row_id = integer(0)
    )
  }

  if (return_best) {
    # Return vector with same length as input, NA for unmatched
    return(dplyr::pull(res, .data$enhanced))
  }

  res <- res |> dplyr::select(-all_of("row_id"))
  res
}
