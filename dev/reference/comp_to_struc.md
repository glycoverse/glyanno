# Convert glycan composition to glycan structure

Given glycan compositions, this function matches them to all compatible
glycan structures in the `glydb` database. Generic, concrete, and mixed
residue identities are matched residue by residue.

## Usage

``` r
comp_to_struc(
  comps,
  db = lifecycle::deprecated(),
  return_best = FALSE,
  glycan_type = NULL,
  species = NULL,
  structure_level = "intact",
  mono_type = "concrete",
  mono_range = NULL
)
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

  **\[deprecated\]** Glycan structures to match against. Can be a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector or any structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  Structures with unresolved floating parts or substituents are excluded
  with a warning. Use the filtering arguments instead. This argument
  cannot be combined with them.

- return_best:

  If `TRUE`, only return the highest confidence match for each
  composition. A custom `db` must have a `confidence` attribute. Default
  is `FALSE`.

- glycan_type:

  Glycan type to select from
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html).

- species:

  Species to select, matched without regard to letter case.

- structure_level:

  Structure resolution, `"intact"` (default) or `"topological"`.

- mono_type:

  Monosaccharide resolution, `"concrete"` (default) or `"generic"`.

- mono_range:

  Named list of monosaccharide count ranges; see
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html).

## Value

If `return_best=TRUE`: A
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector with the same length as `comps`. Unmatched compositions are
returned as `NA`. If `return_best=FALSE`: A tibble with the following
columns:

- `composition`: The glycan compositions, as
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  vector.

- `structure`: The possible glycan structures, as
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector.

- `confidence`: The database confidence score for each structure, or
  `NA` when no score is available. Note that one glycan composition can
  have multiple rows in the result, corresponding to different possible
  glycan structures.

## Details

Filter the built-in database with `glycan_type`, `species`,
`structure_level`, `mono_type`, and `mono_range`.

## See also

[`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html)

## Examples

``` r
comp_to_struc("H5N2")
#> Warning: `db` contains 564 structures with unresolved floating parts or substituents.
#> ℹ Those database structures were excluded from matching.
#> # A tibble: 81 × 3
#>    composition     structure                                          confidence
#>    <comp>          <glydb_st>                                              <dbl>
#>  1 Hex(5)HexNAc(2) Man(b1-2)Man(b1-3)[Man(b1-3)Man(b1-6)]Man(b1-4)Gl…      -1   
#>  2 Hex(5)HexNAc(2) Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)…       5.10
#>  3 Hex(5)HexNAc(2) GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-…      -1   
#>  4 Hex(5)HexNAc(2) Man(a1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)Gl…       3.22
#>  5 Hex(5)HexNAc(2) Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(a1-4)…       1.39
#>  6 Hex(5)HexNAc(2) Man(a1-2)Man(a1-2)Man(a1-3)[Glc(a1-6)]Man(b1-4)Gl…      -1   
#>  7 Hex(5)HexNAc(2) Man(a1-3)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)Gl…      -1   
#>  8 Hex(5)HexNAc(2) Man(a1-3)[Man(a1-6)]Man(a1-6)Man(a1-3)Man(b1-4)Gl…      -1   
#>  9 Hex(5)HexNAc(2) Galf(b1-2)[Man(a1-3)]Man(a1-6)[Man(a1-3)]Man(b1-4…      -1   
#> 10 Hex(5)HexNAc(2) Man(a1-2)Man(a1-3)[Man(a1-6)Man(a1-6)]Man(a1-4)Gl…      -1   
#> # ℹ 71 more rows
```
