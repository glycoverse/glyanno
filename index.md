# glyanno

The goal of glyanno is to enhance the level of information and
resolution for glycans through a hierarchical annotation approach.
Starting from a molecule mass, glyanno calculates possible glycan
compositions. Given a glycan composition, glyanno further deduces
possible glycan structures. For glycan structures with generic
monosaccharides (e.g., “Hex”, “HexNAc”), glyanno refines them into
specific types (e.g., “Glc”, “Gal”). For structures lacking linkage
information (e.g., “Gal(??-?)GalNAc(??-”), glyanno infers the most
likely linkages (e.g., “Gal(b1-3)GalNAc(a1-”).

## Installation

You can install the latest release of glyanno from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("glycoverse/glyanno@*release")
```

Or install the development version:

``` r
remotes::install_github("glycoverse/glyanno")
```

## Documentation

- 🚀 Get started:
  [Here](https://glycoverse.github.io/glyanno/articles/glyanno.html)
- 📚 Reference:
  [Here](https://glycoverse.github.io/glyanno/reference/index.html)

## Role in `glycoverse`

Glyanno enhances the precision of glycan structures, enabling the use of
rough mass spectrometry results for detailed structural analysis. For
instance, the `glyenzy` package requires fully determined glycan
structures. Users can leverage glyanno to make educated guesses from
incomplete data, thereby enabling downstream analysis with glyenzy and
other tools that demand high-resolution structural information.

## Example

``` r
library(glyanno)
library(glyrepr)
library(glydb)

mz_to_comp(406.1325, charge = 1, adduct = "Na+")
#> # A tibble: 4 × 2
#>      mz composition    
#>   <dbl> <comp>         
#> 1  406. Gal(1)GalNAc(1)
#> 2  406. Gal(1)GlcNAc(1)
#> 3  406. Glc(1)GlcNAc(1)
#> 4  406. Man(1)GlcNAc(1)
```

``` r
comp_to_struc("Gal(1)GalNAc(1)")
#> # A tibble: 19 × 2
#>    composition     structure          
#>    <comp>          <struct>           
#>  1 Gal(1)GalNAc(1) Gal(b1-3)GalNAc(a1-
#>  2 Gal(1)GalNAc(1) Gal(b1-3)GalNAc(b1-
#>  3 Gal(1)GalNAc(1) GalNAc(a1-3)Gal(a1-
#>  4 Gal(1)GalNAc(1) Gal(a1-3)GalNAc(b1-
#>  5 Gal(1)GalNAc(1) Gal(b1-4)GalNAc(b1-
#>  6 Gal(1)GalNAc(1) GalNAc(a1-4)Gal(b1-
#>  7 Gal(1)GalNAc(1) GalNAc(a1-2)Gal(b1-
#>  8 Gal(1)GalNAc(1) GalNAc(b1-2)Gal(a1-
#>  9 Gal(1)GalNAc(1) GalNAc(b1-4)Gal(b1-
#> 10 Gal(1)GalNAc(1) Gal(a1-6)GalNAc(a1-
#> 11 Gal(1)GalNAc(1) Gal(b1-6)GalNAc(a1-
#> 12 Gal(1)GalNAc(1) GalNAc(a1-3)Gal(b1-
#> 13 Gal(1)GalNAc(1) GalNAc(b1-3)Gal(b1-
#> 14 Gal(1)GalNAc(1) Gal(b1-6)GalNAc(b1-
#> 15 Gal(1)GalNAc(1) GalNAc(b1-3)Gal(a1-
#> 16 Gal(1)GalNAc(1) GalNAc(b1-4)Gal(a1-
#> 17 Gal(1)GalNAc(1) Gal(a1-3)GalNAc(a1-
#> 18 Gal(1)GalNAc(1) GalNAc(b1-2)Gal(b1-
#> 19 Gal(1)GalNAc(1) Gal(b1-4)GalNAc(a1-
```

``` r
# Use custom db to narrow down the search space
my_db <- glydb_structures(species = "Homo sapiens", glycan_type = "O-GalNAc")
comp_to_struc("Gal(1)GalNAc(1)", db = my_db)
#> # A tibble: 2 × 2
#>   composition     structure          
#>   <comp>          <struct>           
#> 1 Gal(1)GalNAc(1) Gal(b1-3)GalNAc(a1-
#> 2 Gal(1)GalNAc(1) Gal(a1-3)GalNAc(a1-
```
