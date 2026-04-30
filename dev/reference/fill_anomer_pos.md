# Fill anomer positions

Add anomer positions to glycan structures with missing anomer
information based on biological knowledge. For example,
"Gal(??-?)GalNAc(??-" will be converted to "Gal(?1-?)GalNAc(?1-". For
anomer positions that are already specified in the input structures,
this function will not modify them.

## Usage

``` r
fill_anomer_pos(strucs)
```

## Arguments

- strucs:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector of "concrete" monosaccharides.

## Value

A
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector with anomer positions added where missing.

## Details

This function is intended to be used by
[`struc_to_glytoucan()`](https://glycoverse.github.io/glyanno/dev/reference/struc_to_glytoucan.md)
to ensure that glycan structures have complete anomer information before
attempting to inquire the GlyTouCan database.

For these monosaccharides, the anomer position is "2": "Neu5Ac",
"Neu5Gc", "Neu", "Kdn", "Pse", "Leg", "Aci", "4eLeg", "Kdo", "Dha",
"Fru", "Tag", "Sor", and "Psi". All other monosaccharides are assumed to
have anomer positions on "1".

## Examples

``` r
library(glyrepr)
glycans <- as_glycan_structure(c(
  "Gal(??-?)GalNAc(??-",
  "Neu5Ac(??-?)Gal(??-?)GalNAc(??-"
))
fill_anomer_pos(glycans)
#> <glycan_structure[2]>
#> [1] Gal(?1-?)GalNAc(?1-
#> [2] Neu5Ac(?2-?)Gal(?1-?)GalNAc(?1-
#> # Unique structures: 2
```
