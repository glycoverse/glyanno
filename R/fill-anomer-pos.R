#' Fill anomer positions
#'
#' Add anomer positions to glycan structures with missing anomer information based on biological knowledge.
#' For example, "Gal(??-?)GalNAc(??-" will be converted to "Gal(?1-?)GalNAc(?1-".
#' For anomer positions that are already specified in the input structures,
#' this function will not modify them.
#'
#' @details
#' This function is intended to be used by `struc_to_glytoucan()`
#' to ensure that glycan structures have complete anomer information
#' before attempting to inquire the GlyTouCan database.
#'
#' For these monosaccharides, the anomer position is "2":
#' "Neu5Ac", "Neu5Gc", "Neu", "Kdn", "Pse", "Leg", "Aci",
#' "4eLeg", "Kdo", "Dha", "Fru", "Tag", "Sor", and "Psi".
#' All other monosaccharides are assumed to have anomer positions on "1".
#'
#' @param strucs A [glyrepr::glycan_structure()] vector of "concrete" monosaccharides.
#' @returns A [glyrepr::glycan_structure()] vector with anomer positions added where missing.
#'
#' @examples
#' library(glyrepr)
#' glycans <- as_glycan_structure(c(
#'   "Gal(??-?)GalNAc(??-",
#'   "Neu5Ac(??-?)Gal(??-?)GalNAc(??-"
#' ))
#' fill_anomer_pos(glycans)
#'
#' @export
fill_anomer_pos <- function(strucs) {
  checkmate::assert_class(strucs, "glyrepr_structure")
  .assert_concrete(strucs)
  glyrepr::smap_structure(strucs, fill_anomer_pos_single)
}

#' Fill one glycan graph's anomer positions
#'
#' @param struc An igraph glycan structure.
#' @returns An igraph glycan structure with missing anomer positions filled.
#' @noRd
fill_anomer_pos_single <- function(struc) {
  root <- which(igraph::degree(struc, mode = "in") == 0)
  root_mono <- igraph::vertex_attr(struc, "mono", index = root)
  root_anomer <- igraph::graph_attr(struc, "anomer")
  struc <- igraph::set_graph_attr(
    struc,
    "anomer",
    value = fill_anomer_pos_value(root_anomer, root_mono)
  )

  linkages <- igraph::edge_attr(struc, "linkage")
  if (length(linkages) == 0) {
    return(struc)
  }

  edges <- igraph::ends(struc, igraph::E(struc), names = FALSE)
  donor_monos <- igraph::vertex_attr(struc, "mono", index = edges[, 2])
  linkages <- purrr::map2_chr(linkages, donor_monos, fill_anomer_pos_value)
  igraph::set_edge_attr(struc, "linkage", value = linkages)
}

#' Fill the anomer position in one linkage or reducing-end anomer
#'
#' @param anomer A linkage string like `"??-?"` or reducing-end anomer like `"??"`.
#' @param mono A monosaccharide name.
#' @returns `anomer` with the second character filled when missing.
#' @noRd
fill_anomer_pos_value <- function(anomer, mono) {
  if (stringr::str_sub(anomer, 2, 2) != "?") {
    return(anomer)
  }

  stringr::str_c(
    stringr::str_sub(anomer, 1, 1),
    decide_anomer_pos(mono),
    stringr::str_sub(anomer, 3)
  )
}

#' Decide a monosaccharide's anomer position
#'
#' @param mono A monosaccharide name.
#' @returns `"2"` for monosaccharides whose anomer is on carbon 2, otherwise `"1"`.
#' @noRd
decide_anomer_pos <- function(mono) {
  anomer_on_pos2 <- c(
    "Neu5Ac",
    "Neu5Gc",
    "Neu",
    "Kdn",
    "Pse",
    "Leg",
    "Aci",
    "4eLeg",
    "Kdo",
    "Dha",
    "Fru",
    "Tag",
    "Sor",
    "Psi"
  )
  dplyr::if_else(mono %in% anomer_on_pos2, "2", "1")
}
