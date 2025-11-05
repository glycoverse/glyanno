#' Calculate PPM
#'
#' This helper function should be used in [mz_to_comp()] as the `tol` argument.
#' It makes the function using PPM as dynamic tolerance.
#'
#' @param x A numeric scalar of the PPM value.
#'
#' @returns A function that takes a numeric vector of m/z values and returns a numeric vector of tolerances.
#' @examples
#' ppm(10)(2368.84)
#' ppm(10)(c(2368.84, 2368.85))
#' @export
ppm <- function(x) {
  f <- function(mz) x * mz / 1e6
  structure(f, x = x, class = "ppm")
}

#' @export
print.ppm <- function(x, ...) {
  print(glue::glue("ppm({attr(x, 'x')})"))
}