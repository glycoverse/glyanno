# Convert m/z values to glycan composition

Given m/z values, this function matches them to all possible glycan
compositions in the `glydb` database.

## Usage

``` r
mz_to_comp(
  mz,
  tol = ppm(10),
  db = NULL,
  return_best = FALSE,
  charge = 1,
  adduct = "H+",
  mass_dict = NULL
)
```

## Arguments

- mz:

  A numeric vector of m/z values.

- tol:

  A numeric scalar of the tolerance for the m/z value in Da or a
  [`ppm()`](https://glycoverse.github.io/glyanno/dev/reference/ppm.md)
  object for dynamic tolerance. Default is `ppm(10)`.

- db:

  Glycan compositions to match against. Can be a
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  vector or glycan composition strings in Byonic style (e.g.
  Hex(5)HexNAc(2)) or simple style (e.g. H5N4F1S1). If not provided,
  `glydb::glydb_compositions(mono_type = "concrete")` will be used.

- return_best:

  A logical scalar. If `TRUE`, only the match with the highest
  confidence score is returned for each m/z value. The `db` must have a
  `confidence` attribute. Use
  [`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html)
  for `db` to enable this feature. Default is `FALSE`.

- charge:

  Charge to use. Can be 0, 1, 2, 3, -1, -2, -3, etc. 0 means neutral.
  Default is 1.

- adduct:

  Adduct to use. Can be "H+", "K+", "Na+", "NH4+", "H-", "Cl-", "HCO3-".
  Default is "H+".

  - When `charge` is 0, `adduct` is ignored.

  - When `charge` is positive, `adduct` can only be "H+", "K+", "Na+",
    "NH4+".

  - When `charge` is negative, `adduct` can only be "H-", "Cl-",
    "HCO3-".

- mass_dict:

  A named numeric vector of the mass of each monosaccharide residue.
  Default is `glyanno_mass_dict(deriv = "none", mass_type = "mono")`. If
  a custom mass dictionary is provided, please make sure the names of
  the vector are the same as the names in
  [`glyanno_mass_dict()`](https://glycoverse.github.io/glyanno/dev/reference/glyanno_mass_dict.md).

## Value

If `return_best=TRUE`: An unnamed
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
vector with the same length as `mz`. Unmatched m/z values are returned
as `NA`. If `return_best=FALSE`: A tibble with the following columns:

- `mz`: The molecule m/z values, same as the input `mz`.

- `composition`: The possible glycan compositions, as
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  vector.

- `confidence`: The database confidence score for each composition, or
  `NA` when no score is available. Note that one m/z value can have
  multiple rows in the result, corresponding to different possible
  glycan compositions.

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

[`ppm()`](https://glycoverse.github.io/glyanno/dev/reference/ppm.md),
[`glyanno_mass_dict()`](https://glycoverse.github.io/glyanno/dev/reference/glyanno_mass_dict.md)

## Examples

``` r
mz_to_comp(933.3175, charge = 1, adduct = "Na+")
#> # A tibble: 9 × 3
#>      mz composition                    confidence
#>   <dbl> <comp>                              <dbl>
#> 1  933. Glc(1)Gal(2)GlcNAc(1)GalNAc(1)       1.79
#> 2  933. Gal(3)GlcNAc(2)                      1.61
#> 3  933. Glc(1)Gal(2)GlcNAc(2)                1.39
#> 4  933. Man(1)Gal(2)GlcNAc(2)                0   
#> 5  933. Man(3)GlcNAc(2)                      5.61
#> 6  933. Gal(3)GlcNAc(1)GalNAc(1)             1.10
#> 7  933. Glc(1)Gal(2)GalNAc(2)                2.40
#> 8  933. Man(2)Gal(1)GlcNAc(2)               -1   
#> 9  933. Gal(3)GalNAc(2)                     -1   
```
