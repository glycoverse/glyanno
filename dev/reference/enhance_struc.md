# Enhance glycan structure

Given a glycan structure vector of any resolution level (see
[`glyrepr::get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.html)
for details), this function gives all possible glycan structures of
higher resolution level.

## Usage

``` r
enhance_struc(strucs, db = NULL, return_best = FALSE)
```

## Arguments

- strucs:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

- db:

  A
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  The default is
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  at "intact" level. If `db` has a lower or equal resolution level than
  `strucs`, the result will be the same as `strucs` (no enhancement).

- return_best:

  Logical. If `TRUE`, only return the best matching structure (highest
  confidence) for each input structure. `db` must have a `confidence`
  attribute. Default is `FALSE`.

## Value

If `return_best=TRUE`: An unnamed
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector with the same length as `strucs`. Unmatched structures are
returned as `NA`. If `return_best=FALSE`: A tibble with the following
columns:

- `raw`: The original glycan structures.

- `enhanced`: The enhanced glycan structures.

- `confidence`: The database confidence score for each enhanced
  structure, or `NA` when no score is available. Note that one `raw`
  glycan structure can have different `enhanced` glycan structures as
  multiple rows in the result.

## Details

The target resolution level is determined from `db`, defaulting to
[`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
at "intact" level.

## Examples

``` r
# From topological level to intact level
db_intact <- c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
enhance_struc("Gal(??-?)GalNAc(??-", db = db_intact)
#> # A tibble: 2 × 3
#>   raw                 enhanced            confidence
#>   <struct>            <struct>                 <dbl>
#> 1 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(a1-         NA
#> 2 Gal(??-?)GalNAc(??- Gal(b1-4)GalNAc(a1-         NA

# From basic level to topological level
db_topo <- "Gal(??-?)GalNAc(??-"
enhance_struc("Hex(??-?)HexNAc(??-", db = db_topo)
#> # A tibble: 1 × 3
#>   raw                 enhanced            confidence
#>   <struct>            <struct>                 <dbl>
#> 1 Hex(??-?)HexNAc(??- Gal(??-?)GalNAc(??-         NA

# From partial level to intact level
enhance_struc("Gal(b1-?)GalNAc(a1-", db = db_intact)
#> # A tibble: 2 × 3
#>   raw                 enhanced            confidence
#>   <struct>            <struct>                 <dbl>
#> 1 Gal(b1-?)GalNAc(a1- Gal(b1-3)GalNAc(a1-         NA
#> 2 Gal(b1-?)GalNAc(a1- Gal(b1-4)GalNAc(a1-         NA
```
