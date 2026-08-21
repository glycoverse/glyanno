# Enhance glycan structure

Given a glycan structure vector of any resolution level (see
[`glyrepr::get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.html)
for details), this function gives compatible structures with more
specific residue identities or linkage information.

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
  Inputs with unresolved floating parts or substituents are excluded
  with a warning.

- db:

  A
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  Structures with unresolved floating parts or substituents are excluded
  with a warning. The default is
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  at "intact" level.

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
  glycan structure can have different `enhanced` structures as multiple
  rows in the result.

## Details

Input and database vectors may mix generic, concrete, and mixed
residues, as well as topological, partial, and intact structures. Each
database candidate is matched against each input independently.

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

# Refine generic residues without changing the structure level
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
