# Enhance glycan structure

Given a glycan structure vector of any resolution level (see
[`glyrepr::get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.html)
for details), this function gives all possible glycan structures of
higher resolution level.

## Usage

``` r
enhance_struc(strucs, db = NULL, return_best = FALSE, method = "db")
```

## Arguments

- strucs:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

- db:

  A
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  With `method = "db"`, the default is
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  at "intact" level. With `method = "db"`, if `db` has a lower or equal
  resolution level than `strucs`, the result will be the same as
  `strucs` (no enhancement). With `method = "denovo"`, the default is
  [`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
  at "topological" level. A provided `db` cannot be at "basic" level,
  and "partial" or "intact" structures are reduced to "topological"
  before fallback matching.

- return_best:

  Logical. If `TRUE`, only return the best matching structure (highest
  confidence) for each input structure. With `method = "db"`, `db` must
  have a `confidence` attribute. With `method = "denovo"`, this is
  always forced to `TRUE`; an informative message is emitted when
  `FALSE` is supplied. Fallback database candidates require a
  `confidence` attribute. Default is `FALSE`.

- method:

  **\[experimental\]** Enhancement method. `"db"` matches complete
  structures against `db`. `"denovo"` reconstructs topological N-glycans
  from their core and branches.

## Value

With `method = "denovo"`, or if `return_best=TRUE`: An unnamed
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector with the same length as `strucs`. Unmatched structures are
returned as `NA`. With `method = "db"` and `return_best=FALSE`: A tibble
with the following columns:

- `raw`: The original glycan structures.

- `enhanced`: The enhanced glycan structures. Note that one `raw` glycan
  structure can have different `enhanced` glycan structures as multiple
  rows in the result.

## Details

With `method = "db"`, the target resolution level is determined from
`db`, defaulting to
[`glydb::glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
at "intact" level. With `method = "denovo"`, N-glycans are reconstructed
from their core and branches. Inputs that cannot be reconstructed de
novo are matched against a fallback database at "topological" level.

Topological de-novo enhancement preserves optional core fucose and
bisecting GlcNAc residues. Complex N-glycan branches are matched against
the internal branch data. For hybrid N-glycans, the all-Hex arm is
assigned as mannose and the HexNAc-bearing arm uses branch matching.
High-mannose candidates are constrained to core-aligned subtrees of the
Glc3Man9 precursor.

## Examples

``` r
# From topological level to intact level
db_intact <- c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-")
enhance_struc("Gal(??-?)GalNAc(??-", db = db_intact)
#> # A tibble: 2 × 2
#>   raw                 enhanced           
#>   <struct>            <struct>           
#> 1 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(a1-
#> 2 Gal(??-?)GalNAc(??- Gal(b1-4)GalNAc(a1-

# From basic level to topological level
db_topo <- "Gal(??-?)GalNAc(??-"
enhance_struc("Hex(??-?)HexNAc(??-", db = db_topo)
#> # A tibble: 1 × 2
#>   raw                 enhanced           
#>   <struct>            <struct>           
#> 1 Hex(??-?)HexNAc(??- Gal(??-?)GalNAc(??-

# From partial level to intact level
enhance_struc("Gal(b1-?)GalNAc(a1-", db = db_intact)
#> # A tibble: 2 × 2
#>   raw                 enhanced           
#>   <struct>            <struct>           
#> 1 Gal(b1-?)GalNAc(a1- Gal(b1-3)GalNAc(a1-
#> 2 Gal(b1-?)GalNAc(a1- Gal(b1-4)GalNAc(a1-

# De-novo enhancement of an N-glycan
n_basic <- glyrepr::n_glycan_core(linkage = FALSE, mono_type = "generic")
enhance_struc(
  n_basic,
  method = "denovo",
  return_best = TRUE
)
#> <glycan_structure[1]>
#> [1] Man(??-?)[Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-
#> # Unique structures: 1
```
