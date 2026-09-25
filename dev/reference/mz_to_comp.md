# Convert m/z values to glycan composition

Given m/z values, this function matches them to all possible glycan
compositions in the `glydb` database.

## Usage

``` r
mz_to_comp(
  mz,
  tol = ppm(10),
  db = lifecycle::deprecated(),
  return_best = FALSE,
  charge = 1,
  adduct = "H+",
  mass_dict = NULL,
  glycan_type = NULL,
  species = NULL,
  mono_type = "concrete",
  mono_range = NULL
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

  **\[deprecated\]** Glycan compositions to match against. Can be a
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  vector or glycan composition strings in Byonic style (e.g.
  Hex(5)HexNAc(2)) or simple style (e.g. H5N4F1S1). Use the filtering
  arguments instead. This argument cannot be combined with them.

- return_best:

  A logical scalar. If `TRUE`, only the match with the highest
  confidence score is returned for each m/z value. A custom `db` must
  have a `confidence` attribute. Default is `FALSE`.

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

- glycan_type:

  Glycan type to select from
  [`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html).

- species:

  Species to select, matched without regard to letter case.

- mono_type:

  Monosaccharide resolution, `"concrete"` (default) or `"generic"`.

- mono_range:

  Named list of monosaccharide count ranges; see
  [`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html).

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

## Filter the built-in database

Use `glycan_type`, `species`, `mono_type`, and `mono_range` to select
compositions from
[`glydb::glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html).
For example,

    mz_to_comp(mz, species = "Homo sapiens", glycan_type = "N")

## See also

[`ppm()`](https://glycoverse.github.io/glyanno/dev/reference/ppm.md),
[`glyanno_mass_dict()`](https://glycoverse.github.io/glyanno/dev/reference/glyanno_mass_dict.md)

## Examples

``` r
mz_to_comp(933.3175, charge = 1, adduct = "Na+")
#> # A tibble: 10 × 3
#>       mz composition                    confidence
#>    <dbl> <glydb_cm>                          <dbl>
#>  1  933. Man(3)GlcNAc(2)                      5.61
#>  2  933. Glc(1)Gal(2)GlcNAc(1)GalNAc(1)       1.79
#>  3  933. Gal(3)GlcNAc(2)                      1.61
#>  4  933. Glc(1)Gal(2)GlcNAc(2)                1.39
#>  5  933. Man(1)Gal(2)GlcNAc(2)                0   
#>  6  933. Gal(3)GlcNAc(1)GalNAc(1)             1.10
#>  7  933. Glc(1)Gal(2)GalNAc(2)                2.40
#>  8  933. Man(2)Gal(1)GlcNAc(2)                1.39
#>  9  933. Glc(1)GlcNAc(2)Galf(2)              -1   
#> 10  933. Gal(3)GalNAc(2)                     -1   
```
