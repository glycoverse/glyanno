# Reconstruct topological N-glycan structures de novo

**\[experimental\]**

Reconstructs basic-resolution N-glycans from their core and branches.
Inputs that cannot be reconstructed de novo are matched against a
fallback database at "topological" level.

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

- fallback_db:

  A
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  The default is
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  at "topological" level. A provided database cannot be at "basic"
  level. "partial" or "intact" structures are reduced to "topological"
  before fallback matching. Fallback candidates require a `confidence`
  attribute.

## Value

An unnamed
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector with the same length as `strucs`. Unmatched structures are
returned as `NA`.

## Examples

``` r
n_basic <- glyrepr::n_glycan_core(linkage = FALSE, mono_type = "generic")
enhance_struc_denovo(n_basic)
#> <glycan_structure[1]>
#> [1] Man(??-?)[Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-
#> # Unique structures: 1
```
