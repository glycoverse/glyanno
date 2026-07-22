# Assign GlyTouCan accessions to glycan structures

This function takes a vector of glycan structures and returns the
corresponding GlyTouCan accessions. It first looks up structures from
[glydb::glydb_data](https://glycoverse.github.io/glydb/reference/glydb_data.html),
then uses the GlycanFormatConverter API maintained by the Glycosmos
project for structures not found locally.

## Usage

``` r
struc_to_glytoucan(strucs)
```

## Arguments

- strucs:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of glycan text representations supported
  by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  The glycan structure must have "concrete" monosaccharides (e.g., Gal,
  GalNAc).

## Value

A character vector of GlyTouCan accessions corresponding to the input
glycan structures. If a structure cannot be converted to a GlyTouCan
accession, the corresponding entry will be `NA`.

## Details

For "topological" structures (e.g., "Gal(??-?)GalNAc(??-"), this
function will first call
[`glyrepr::fill_anomer_pos()`](https://glycoverse.github.io/glyrepr/reference/fill_anomer_pos.html)
to fill in the missing anomeric positions before looking up accessions.
This is necessary because all glycan structures in GlyTouCan must have
defined anomeric positions.
