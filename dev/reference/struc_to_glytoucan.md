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

## Value

A character vector of GlyTouCan accessions corresponding to the input
glycan structures. If a structure cannot be converted to a GlyTouCan
accession, the corresponding entry will be `NA`.
