test_that("Neu and Kdn dictionary entries match generic residue names", {
  compositions <- glyrepr::as_glycan_composition(c(
    "Neu(1)",
    "Kdn(1)",
    "Neu(1)Kdn(1)"
  ))
  generic <- glyrepr::convert_to_generic(compositions)

  for (deriv in c("none", "permethyl", "peracetyl")) {
    for (mass_type in c("mono", "average")) {
      dictionary <- glyanno_mass_dict(deriv, mass_type)
      expected <- unname(c(
        dictionary["Neu"],
        dictionary["Kdn"],
        dictionary["Neu"] + dictionary["Kdn"]
      ) + dictionary["red_end"])
      expect_equal(
        calculate_mz(compositions, charge = 0, mass_dict = dictionary),
        expected
      )
      expect_equal(
        calculate_mz(generic, charge = 0, mass_dict = dictionary),
        expected
      )
      db <- list(
        ids = c(3L, 1L, 2L),
        view = list(masses = list(), mass_counts = .mass_count_matrix(generic))
      )
      expect_equal(.cached_mz(db, 0, "H+", dictionary), expected[db$ids])
      expect_false(any(c("gNeu", "gKdn") %in% names(dictionary)))
    }
  }
})

test_that("generic aliases retain custom Neu and Kdn mass values", {
  dictionary <- glyanno_mass_dict()
  dictionary[c("Neu", "Kdn")] <- c(300, 400)
  expect_equal(
    calculate_mz("gNeu(1)gKdn(1)", charge = 0, mass_dict = dictionary),
    unname(700 + dictionary["red_end"])
  )
})
