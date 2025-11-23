# Enhance glycan composition

Given a generic glycan composition (e.g. Hex(5)HexNAc(2)), this function
gives all possible concrete glycan compositions (e.g. Man(5)GlcNAc(2)).

## Usage

``` r
enhance_comp(comps, db = NULL)
```

## Arguments

- comps:

  A
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  vector, or a character vector of glycan composition strings of Byonic
  or simple style (e.g. "Hex(5)HexNAc(2)", "H5N4F1S1"). Generic
  compositions (e.g. Hex(5)HexNAc(2)) will be matched to all possible
  concrete compositions in `db`. Concrete compositions (e.g.
  Man(5)GlcNAc(2)) will be returned as is.

- db:

  A
  [`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html)
  vector, or a character vector of glycan composition strings of Byonic
  or simple style (e.g. "Man(5)GlcNAc(2)", "H5N4F1S1"). All compositions
  in `db` must be concrete (e.g. Man(5)GlcNAc(2)). If not provided,
  `glydb::glydb_compositions(mono_type = "concrete")` will be used.

## Value

A tibble with the following columns:

- `raw`: The original compositions.

- `enhanced`: The enhanced compositions. Note that one `raw` composition
  can have different `enhanced` compositions as multiple rows in the
  result.

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
enhance_comp("Hex(5)HexNAc(2)")
#> # A tibble: 5 × 2
#>   raw             enhanced             
#>   <comp>          <comp>               
#> 1 Hex(5)HexNAc(2) Man(5)GlcNAc(2)      
#> 2 Hex(5)HexNAc(2) Glc(1)Man(4)GlcNAc(2)
#> 3 Hex(5)HexNAc(2) Man(4)Gal(1)GlcNAc(2)
#> 4 Hex(5)HexNAc(2) Glc(1)Gal(4)GlcNAc(2)
#> 5 Hex(5)HexNAc(2) Man(3)Gal(2)GlcNAc(2)
```
