# Enhance glycan structure

Given a glycan structure vector of any resolution level (see
[`glyrepr::get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.html)
for details), this function gives compatible structures with more
specific residue identities or linkage information.

## Usage

``` r
enhance_struc(
  strucs,
  db = lifecycle::deprecated(),
  return_best = FALSE,
  glycan_type = NULL,
  species = NULL,
  structure_level = "intact",
  mono_type = "concrete",
  mono_range = NULL
)
```

## Arguments

- strucs:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  Inputs with unresolved floating parts or substituents are excluded
  with a warning.

- db:

  **\[deprecated\]** A
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  Structures with unresolved floating parts or substituents are excluded
  with a warning. Use the filtering arguments instead. This argument
  cannot be combined with them.

- return_best:

  Logical. If `TRUE`, only return the best matching structure (highest
  confidence) for each input structure. A custom `db` must have a
  `confidence` attribute. Default is `FALSE`.

- glycan_type:

  Glycan type to select from
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html).

- species:

  Species to select, matched without regard to letter case.

- structure_level:

  Structure resolution, `"intact"` (default) or `"topological"`.

- mono_type:

  Monosaccharide resolution, `"concrete"` (default) or `"generic"`.

- mono_range:

  Named list of monosaccharide count ranges; see
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html).

## Value

If `return_best=TRUE`: An unnamed
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector with the same length as `strucs`. Unmatched structures are
returned as `NA`. If `return_best=FALSE`: A tibble with the following
columns:

- `raw`: The original glycan structures.

- `enhanced`: The enhanced glycan structures.

- `confidence`: The database confidence score for each enhanced
  structure, or `NA` when no score is available. Note that one `raw`
  glycan structure can have different `enhanced` structures as multiple
  rows in the result.

## Details

Input and database vectors may mix generic, concrete, and mixed
residues, as well as topological, partial, and intact structures. Each
database candidate is matched against each input independently.

## Examples

``` r
enhance_struc("Gal(??-?)GalNAc(??-", glycan_type = "O-GalNAc",
  species = "Homo sapiens")
#> Warning: `db` contains 9 structures with unresolved floating parts or substituents.
#> ℹ Those database structures were excluded from matching.
#> # A tibble: 3 × 3
#>   raw                 enhanced            confidence
#>   <struct>            <glydb_st>               <dbl>
#> 1 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(a1-       5.32
#> 2 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(b1-       1.10
#> 3 Gal(??-?)GalNAc(??- Gal(a1-3)GalNAc(a1-       1.61
```
