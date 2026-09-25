# Enhance glycan composition

Given a generic or mixed glycan composition (e.g. Hex(5)HexNAc(2)), this
function gives all possible compatible concrete glycan compositions
(e.g. Man(5)GlcNAc(2)).

## Usage

``` r
enhance_comp(
  comps,
  db = lifecycle::deprecated(),
  return_best = FALSE,
  glycan_type = NULL,
  species = NULL,
  mono_type = "concrete",
  mono_range = NULL
)
```

## Arguments

- comps:

  A
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  vector, or a character vector of glycan composition strings of Byonic
  or simple style (e.g. "Hex(5)HexNAc(2)", "H5N4F1S1"). Generic and
  mixed compositions are matched to compatible concrete compositions in
  `db`. Concrete compositions are returned as is.

- db:

  **\[deprecated\]** A
  [`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html)
  vector, or a character vector of glycan composition strings of Byonic
  or simple style (e.g. "Man(5)GlcNAc(2)", "H5N4F1S1"). All compositions
  in `db` must be concrete (e.g. Man(5)GlcNAc(2)). Use the filtering
  arguments instead. This argument cannot be combined with them.

- return_best:

  Logical. If `TRUE`, only return the highest confidence match for each
  input composition. A custom `db` must have a `confidence` attribute.
  Defaults to `FALSE`.

- glycan_type:

  Glycan type to select from
  [`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html).

- species:

  Species to select, matched without regard to letter case.

- mono_type:

  Monosaccharide resolution; only `"concrete"` is supported for
  enhancement.

- mono_range:

  Named list of monosaccharide count ranges; see
  [`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html).

## Value

If `return_best=TRUE`: A
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
vector with the same length as `comps`. Unmatched compositions are
returned as `NA`. If `return_best=FALSE`: A tibble with the following
columns:

- `raw`: The original compositions.

- `enhanced`: The enhanced compositions.

- `confidence`: The database confidence score for each enhanced
  composition, or `NA` when no score is available. Note that one `raw`
  composition can have different `enhanced` compositions as multiple
  rows in the result.

## Details

Filter the built-in database with `glycan_type`, `species`, `mono_type`,
and `mono_range`.

## Examples

``` r
enhance_comp("Hex(5)HexNAc(2)")
#> # A tibble: 6 × 3
#>   raw             enhanced               confidence
#>   <comp>          <glydb_cm>                  <dbl>
#> 1 Hex(5)HexNAc(2) Man(5)GlcNAc(2)             5.10 
#> 2 Hex(5)HexNAc(2) Glc(1)Gal(4)GlcNAc(2)       1.10 
#> 3 Hex(5)HexNAc(2) Man(4)GlcNAc(2)Galf(1)     -1    
#> 4 Hex(5)HexNAc(2) Man(4)Gal(1)GlcNAc(2)       0    
#> 5 Hex(5)HexNAc(2) Man(3)Gal(2)GlcNAc(2)       0.693
#> 6 Hex(5)HexNAc(2) Glc(1)Man(4)GlcNAc(2)       0.693
```
