# Assign GlyTouCan accessions to glycan structures

This function takes a vector of glycan structures and returns the
corresponding GlyTouCan accessions. Under the hood, it uses the
GlycanFormatConverter API maintained by the Glycosmos project.

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
[`fill_anomer_pos()`](https://glycoverse.github.io/glyanno/dev/reference/fill_anomer_pos.md)
to fill in the missing anomeric positions before querying the API. This
is necessary because all glycan structures in GlyTouCan must have
defined anomeric positions.
