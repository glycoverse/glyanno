test_that("de-novo enhancement preserves substituents on core additions", {
  withr::local_options(lifecycle_verbosity = "quiet")
  cases <- c(
    paste0(
      "HexNAc(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]",
      "[HexNAc6S(??-?)]Hex(??-?)HexNAc(??-?)[dHex3S(??-?)]HexNAc(??-"
    ),
    paste0(
      "Hex(??-?)Hex(??-?)[Hex(??-?)][HexNAc6S(??-?)]Hex(??-?)",
      "HexNAc(??-?)[dHex3S(??-?)]HexNAc(??-"
    )
  )
  inputs <- glyrepr::as_glycan_structure(cases)
  result <- enhance_struc_denovo(inputs)

  expect_false(anyNA(result))
  expect_equal(
    glyrepr::as_glycan_composition(glyrepr::convert_to_generic(result)),
    glyrepr::as_glycan_composition(inputs)
  )
  for (i in seq_along(result)) {
    nodes <- glyrepr::structure_nodes(result[i])
    expect_equal(nodes$sub[nodes$mono == "Fuc"], "3S")
    expect_true(any(nodes$mono == "GlcNAc" & nodes$sub == "6S"))
  }
})
