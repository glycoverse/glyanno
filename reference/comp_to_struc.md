# Convert glycan composition to glycan structure

Given glycan compositions, this function matches them to all possible
glycan structures in the `glydb` database.

## Usage

``` r
comp_to_struc(comps, db = NULL, return_best = FALSE)
```

## Arguments

- comps:

  Glycan compositions to match against. Can be either:

  - A
    [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
    vector.

  - Byonic style composition strings (e.g. Hex(5)HexNAc(2)).

  - Simple style composition strings (e.g. H5N4F1S1).

- db:

  Glycan structures to match against. Can be a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector or any structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  If not provided, `glydb::glydb_structures(structure_level = "intact")`
  will be used.

- return_best:

  If `TRUE`, only return the highest confidence match for each
  composition. Requires `db` to have a `confidence` attribute. Use
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  for `db` to enable this feature. Default is `FALSE`.

## Value

A tibble with the following columns:

- `composition`: The glycan compositions, as
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  vector.

- `structure`: The possible glycan structures, as
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector. Note that one glycan composition can have multiple rows in the
  result, corresponding to different possible glycan structures.

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

## See also

[`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html)

## Examples

``` r
comp_to_struc("H5N2")
#> # A tibble: 79 × 2
#>    composition     structure                                                    
#>    <comp>          <struct>                                                     
#>  1 Hex(5)HexNAc(2) Man(b1-2)Man(b1-3)[Man(b1-3)Man(b1-6)]Man(b1-4)GlcNAc(b1-4)G…
#>  2 Hex(5)HexNAc(2) GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4…
#>  3 Hex(5)HexNAc(2) Man(a1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)G…
#>  4 Hex(5)HexNAc(2) Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(a1-4)GlcNAc(b1-4…
#>  5 Hex(5)HexNAc(2) Man(a1-2)Man(a1-2)Man(a1-3)[Glc(a1-6)]Man(b1-4)GlcNAc(b1-4)G…
#>  6 Hex(5)HexNAc(2) Man(a1-3)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)G…
#>  7 Hex(5)HexNAc(2) Man(a1-3)[Man(a1-6)]Man(a1-6)Man(a1-3)Man(b1-4)GlcNAc(b1-4)G…
#>  8 Hex(5)HexNAc(2) Man(a1-2)Man(a1-3)[Man(a1-6)Man(a1-6)]Man(a1-4)GlcNAc(b1-4)G…
#>  9 Hex(5)HexNAc(2) Man(a1-2)Man(a1-3)Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)G…
#> 10 Hex(5)HexNAc(2) Gal(a1-6)Man(b1-3)[Man(b1-6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)G…
#> # ℹ 69 more rows
```
