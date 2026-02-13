test_that(".check_confidence_attr passes when return_best is FALSE", {
  db <- glyrepr::glycan_composition(c(Hex = 1))
  expect_no_error(glyanno:::.check_confidence_attr(db, FALSE))
})

test_that(".check_confidence_attr passes when return_best is TRUE and confidence exists", {
  db <- glyrepr::glycan_composition(c(Hex = 1))
  attr(db, "confidence") <- 1.0
  expect_no_error(glyanno:::.check_confidence_attr(db, TRUE))
})

test_that(".check_confidence_attr errors when return_best is TRUE and no confidence", {
  db <- glyrepr::glycan_composition(c(Hex = 1))
  expect_error(
    glyanno:::.check_confidence_attr(db, TRUE),
    "must have a .*confidence.* attribute"
  )
})
