#' Assign GlyTouCan accessions to glycan structures
#'
#' This function takes a vector of glycan structures and returns the corresponding GlyTouCan accessions.
#' Under the hood, it uses the GlycanFormatConverter API maintained by the Glycosmos project.
#'
#' @details
#' For "topological" structures (e.g., "Gal(??-?)GalNAc(??-"),
#' this function will first call [fill_anomer_pos()] to fill in the missing anomeric positions before querying the API.
#' This is necessary because all glycan structures in GlyTouCan must have defined anomeric positions.
#'
#' @param strucs A [glyrepr::glycan_structure()] vector,
#'   or a character vector of glycan text representations supported by [glyparse::auto_parse()].
#'
#' @returns A character vector of GlyTouCan accessions corresponding to the input glycan structures.
#'   If a structure cannot be converted to a GlyTouCan accession, the corresponding entry will be `NA`.
#' @export
struc_to_glytoucan <- function(strucs) {
  strucs <- .ensure_glycan_structure(strucs)
  strucs <- fill_anomer_pos(strucs)
  missing_strucs <- is.na(strucs)
  accessions <- rep(NA_character_, length(strucs))

  if (all(missing_strucs)) {
    return(accessions)
  }

  iupacs <- glyrepr::structure_to_iupac(strucs[!missing_strucs])
  base_url <- "https://api.glycosmos.org/glycanformatconverter/2.8.2/iupaccondensed2wurcs/"
  encoded_iupacs <- purrr::map_chr(
    iupacs,
    utils::URLencode,
    reserved = TRUE,
    repeated = TRUE
  )
  urls <- paste0(base_url, encoded_iupacs)

  reqs <- purrr::map(urls, httr2::request) |>
    purrr::map(httr2::req_throttle, capacity = 60) |>
    purrr::map(httr2::req_retry, max_tries = 2)
  resps <- purrr::map(reqs, httr2::req_perform)

  accessions[!missing_strucs] <- purrr::map_chr(resps, function(resp) {
    if (httr2::resp_status(resp) == 200) {
      content <- httr2::resp_body_json(resp)
      if (!is.null(content$id)) {
        return(content$id)
      }
    }
    NA_character_
  })
  accessions
}
