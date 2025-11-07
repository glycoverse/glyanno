test_that("ppm works", {
  tol <- ppm(10)
  expect_snapshot(tol)
  expect_equal(tol(2368.84), 0.0236884)
})

test_that("ppm works for vector input", {
  expect_equal(ppm(10)(c(1000, 2000)), c(0.01, 0.02))
})