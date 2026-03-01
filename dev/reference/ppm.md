# Calculate PPM

This helper function should be used in
[`mz_to_comp()`](https://glycoverse.github.io/glyanno/dev/reference/mz_to_comp.md)
as the `tol` argument. It makes the function using PPM as dynamic
tolerance.

## Usage

``` r
ppm(x)
```

## Arguments

- x:

  A numeric scalar of the PPM value.

## Value

A function that takes a numeric vector of m/z values and returns a
numeric vector of tolerances.

## Examples

``` r
ppm(10)(2368.84)
#> [1] 0.0236884
ppm(10)(c(2368.84, 2368.85))
#> [1] 0.0236884 0.0236885
```
