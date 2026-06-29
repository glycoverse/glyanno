# Calculate m/z values of glycans

This function calculates m/z values from glycans. Different adducts,
mass types, and derivatives are supported. Custom mass dictionaries are
also supported.

## Usage

``` r
calculate_mz(glycans, charge = 1, adduct = "H+", mass_dict = NULL, safe = TRUE)
```

## Arguments

- glycans:

  Glycans to calculate m/z values from. Valid inputs include:

  - A
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    vector.

  - A
    [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
    vector.

  - Byonic style composition strings (e.g. Hex(5)HexNAc(2)).

  - Simple style composition strings (e.g. H5N4F1S1).

  - Any structure strings supported by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

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

- safe:

  Whether to raise an error when unsupported monosaccharides or
  substituents are found.

  - If `TRUE` (default), an error will be raised.

  - If `FALSE`, a warning will be raised and m/z values for invalid
    glycans will be set to NA.

## Value

A numeric vector of m/z values.

## See also

[`glyanno_mass_dict()`](https://glycoverse.github.io/glyanno/dev/reference/glyanno_mass_dict.md)

## Examples

``` r
library(glyrepr)
#> 
#> Attaching package: ‘glyrepr’
#> The following object is masked from ‘package:glyanno’:
#> 
#>     fill_anomer_pos

# Different input types
calculate_mz(glycan_composition(c(Gal = 1, GalNAc = 1)), charge = 0)
#> [1] 383.1428
calculate_mz("Gal(1)GalNAc(1)", charge = 0)
#> [1] 383.1428
calculate_mz(as_glycan_structure("Gal(b1-3)GalNAc(a1-"), charge = 0)
#> [1] 383.1428
calculate_mz("Gal(b1-3)GalNAc(a1-", charge = 0)
#> [1] 383.1428

# For common situation in MALDI-TOF MS
calculate_mz("Man(5)GlcNAc(2)", charge = 1,adduct = "Na+")
#> [1] 1257.423

# For common situation in ESI MS
calculate_mz("Man(5)GlcNAc(2)", charge = 1, adduct = "H+")
#> [1] 1235.441

# Calculate permethylated m/z values
calculate_mz("Man(5)GlcNAc(2)", charge = 0, mass_dict = glyanno_mass_dict(deriv = "permethyl"))
#> [1] 1556.793

# Calculate average mass
calculate_mz("Man(5)GlcNAc(2)", charge = 0, mass_dict = glyanno_mass_dict(mass_type = "average"))
#> [1] 1235.117

# Vectorization
calculate_mz(c("Man(5)GlcNAc(2)", "Gal(1)GalNAc(1)"), charge = 1, adduct = "Na+")
#> [1] 1257.4231  406.1325
```
