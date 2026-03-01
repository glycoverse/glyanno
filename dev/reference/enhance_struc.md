# Enhance glycan structure

Given a glycan structure of any resolution level (see
[`glyrepr::get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.html)
for details), this function gives all possible glycan structures of
higher resolution level.

## Usage

``` r
enhance_struc(strucs, to_level = "intact", db = NULL, return_best = FALSE)
```

## Arguments

- strucs:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  Glycan structures with level higher or same as `to_level` will be
  returned as is. Glycan structures with level lower than `to_level`
  will be enhanced to the level of `to_level`.

- to_level:

  The resolution level to enhance to. Can be "intact" or "topological".
  Default is "intact".

- db:

  A
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  All structures in `db` must be at the same resolution level as
  `to_level`. If not provided, a default structure vector is loaded from
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  with the specified `to_level`.

- return_best:

  Logical. If `TRUE`, only return the best matching structure (highest
  confidence) for each input structure. Requires `db` to have a
  `confidence` attribute. Use
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  for `db` to enable this feature. Default is `FALSE`.

## Value

A tibble with the following columns:

- `raw`: The original glycan structures.

- `enhanced`: The enhanced glycan structures. Note that one `raw` glycan
  structure can have different `enhanced` glycan structures as multiple
  rows in the result.

## How to set `db`

The `db` parameter is very important for all functions in this package.
By default, it uses all available glycans in the `glydb` package, which
is usually larger than what you need. You can use helper functions in
`glydb` to narrow down the database, e.g.
[`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html)
or
[`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html).

You can use the `species` and `glycan_type` parameters to focus on
specific species and glycan type. For example, if you are only
interested in N-glycan compositions in human, you can use
`glydb::glydb_compositions(species = "Homo sapiens", glycan_type = "N")`.
Also, you can decide the level of information in the database by setting
`mono_type` of
[`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html)
and `structure_level` of
[`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html).

You can then pass the result to the `db` parameter of this function. For
example,

    my_db <- glydb::glydb_compositions(species = "Homo sapiens", glycan_type = "N")
    mz_to_comp(mz, db = my_db)

## Examples

``` r
# From topological level to intact level
enhance_struc("Gal(??-?)GalNAc(??-", to_level = "intact")
#> # A tibble: 9 × 2
#>   raw                 enhanced           
#>   <struct>            <struct>           
#> 1 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(a1-
#> 2 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(b1-
#> 3 Gal(??-?)GalNAc(??- Gal(a1-3)GalNAc(b1-
#> 4 Gal(??-?)GalNAc(??- Gal(b1-4)GalNAc(b1-
#> 5 Gal(??-?)GalNAc(??- Gal(a1-6)GalNAc(a1-
#> 6 Gal(??-?)GalNAc(??- Gal(b1-6)GalNAc(a1-
#> 7 Gal(??-?)GalNAc(??- Gal(b1-6)GalNAc(b1-
#> 8 Gal(??-?)GalNAc(??- Gal(a1-3)GalNAc(a1-
#> 9 Gal(??-?)GalNAc(??- Gal(b1-4)GalNAc(a1-

# From basic level to topological level
enhance_struc("Hex(??-?)HexNAc(??-", to_level = "topological")
#> # A tibble: 4 × 2
#>   raw                 enhanced           
#>   <struct>            <struct>           
#> 1 Hex(??-?)HexNAc(??- Gal(??-?)GalNAc(??-
#> 2 Hex(??-?)HexNAc(??- Gal(??-?)GlcNAc(??-
#> 3 Hex(??-?)HexNAc(??- Glc(??-?)GlcNAc(??-
#> 4 Hex(??-?)HexNAc(??- Man(??-?)GlcNAc(??-

# From partial level to intact level
enhance_struc("Gal(b1-?)GalNAc(a1-", to_level = "intact")
#> # A tibble: 3 × 2
#>   raw                 enhanced           
#>   <struct>            <struct>           
#> 1 Gal(b1-?)GalNAc(a1- Gal(b1-3)GalNAc(a1-
#> 2 Gal(b1-?)GalNAc(a1- Gal(b1-6)GalNAc(a1-
#> 3 Gal(b1-?)GalNAc(a1- Gal(b1-4)GalNAc(a1-
```
