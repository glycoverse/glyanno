test_that(".check_confidence_attr passes when return_best is FALSE", {
  expect_no_error(glyanno:::.check_confidence_attr(NULL, FALSE))
})

test_that(".check_confidence_attr passes when return_best is TRUE and confidence exists", {
  confidences <- c(1.0, 2.0)
  expect_no_error(glyanno:::.check_confidence_attr(confidences, TRUE))
})

test_that(".check_confidence_attr errors when return_best is TRUE and no confidence", {
  expect_error(
    glyanno:::.check_confidence_attr(NULL, TRUE),
    "must have a .*confidence.* attribute"
  )
})
