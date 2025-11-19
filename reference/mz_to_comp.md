# Convert m/z values to glycan composition

Given m/z values, this function matches them to all possible glycan
compositions in the `glydb` database.

## Usage

``` r
mz_to_comp(
  mz,
  tol = ppm(10),
  db = NULL,
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
  [`ppm()`](https://glycoverse.github.io/glyanno/reference/ppm.md)
  object for dynamic tolerance. Default is `ppm(10)`.

- db:

  Glycan compositions to match against. Can be a
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  vector or glycan composition strings in Byonic style (e.g.
  Hex(5)HexNAc(2)) or simple style (e.g. H5N4F1S1). If not provided,
  `glydb::glydb_compositions(mono_type = "concrete")` will be used.

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
  [`glyanno_mass_dict()`](https://glycoverse.github.io/glyanno/reference/glyanno_mass_dict.md).

## Value

A tibble with the following columns:

- `mz`: The molecule m/z values, same as the input `mz`.

- `composition`: The possible glycan compositions, as
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
  vector. Note that one m/z value can have multiple rows in the result,
  corresponding to different possible glycan compositions.

## See also

[`ppm()`](https://glycoverse.github.io/glyanno/reference/ppm.md),
[`glyanno_mass_dict()`](https://glycoverse.github.io/glyanno/reference/glyanno_mass_dict.md)

## Examples

``` r
mz_to_comp(933.3175, charge = 1, adduct = "Na+")
#> # A tibble: 9 × 2
#>      mz composition                   
#>   <dbl> <comp>                        
#> 1  933. Glc(1)Gal(2)GlcNAc(1)GalNAc(1)
#> 2  933. Man(3)GlcNAc(2)               
#> 3  933. Gal(3)GlcNAc(2)               
#> 4  933. Gal(3)GlcNAc(1)GalNAc(1)      
#> 5  933. Glc(1)Gal(2)GalNAc(2)         
#> 6  933. Man(1)Gal(2)GlcNAc(2)         
#> 7  933. Glc(1)Gal(2)GlcNAc(2)         
#> 8  933. Man(2)Gal(1)GlcNAc(2)         
#> 9  933. Gal(3)GalNAc(2)               
```
