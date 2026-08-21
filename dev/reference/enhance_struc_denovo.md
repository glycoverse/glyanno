# Reconstruct topological N-glycan structures de novo

**\[experimental\]**

Reconstructs generic topological N-glycans from their core and branches.
Every non-missing input must have only generic residues and no linkage
information. Inputs that cannot be reconstructed de novo are matched
against a fallback database at "topological" level.

Topological de-novo enhancement preserves optional core fucose and
bisecting GlcNAc residues. Complex N-glycan branches are matched against
the internal branch data. For hybrid N-glycans, the all-Hex arm is
assigned as mannose and the HexNAc-bearing arm uses branch matching.
High-mannose candidates are constrained to core-aligned subtrees of the
Glc3Man9 precursor.

## Usage

``` r
enhance_struc_denovo(strucs, fallback_db = NULL)
```

## Arguments

- strucs:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  Every non-missing element must be generic and topological. Missing
  values are allowed and preserved. Inputs with unresolved floating
  parts or substituents are replaced with missing values and produce a
  warning.

- fallback_db:

  A
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  The default is
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  at "topological" level. Every non-missing candidate must have concrete
  residues. Linkage information is removed before fallback matching.
  Candidates with unresolved floating parts or substituents are excluded
  with a warning. Fallback candidates require a `confidence` attribute.

## Value

An unnamed
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector with the same length as `strucs`. Unmatched structures are
returned as `NA`.

## Examples

``` r
n_generic <- glyrepr::n_glycan_core(linkage = FALSE, mono_type = "generic")
enhance_struc_denovo(n_generic)
#> <glycan_structure[1]>
#> [1] Man(??-?)[Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-
#> # Unique structures: 1
```
