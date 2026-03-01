# Get Started with glyanno

## The Challenge

Mass spectrometry-based glycomics and glycoproteomics often yield data
with limited structural resolution. For example, MALDI-TOF MS-based
glycomics typically provides only m/z values, while many LC-MS/MS-based
glycoproteomics workflows identify basic glycan compositions but lack
isomer specificity. Mass spectrometry alone frequently fails to
distinguish specific monosaccharide isomers (e.g., differentiating
Mannose from Galactose), or resolve linkage details (e.g., a1-3 vs
b1-4).

In theory, full structural resolution is achievable by combining MS with
orthogonal techniques such as enzymatic digestion or NMR; however, these
approaches are often resource-intensive and time-consuming.
Unfortunately, advanced downstream analyses—such as the glycan
biosynthetic pathway reconstruction offered by
[glyenzy](https://github.com/glycoverse/glyenzy)— require fully resolved
glycan structures, including specific monosaccharide identities and
linkage information.

## What is glyanno?

`glyanno` bridges this gap. This package is designed to maximize the
utility of your mass spectrometry data by annotating it with probable
biological context derived from knowledge bases. The package features
four core functions:

- [`mz_to_comp()`](https://glycoverse.github.io/glyanno/dev/reference/mz_to_comp.md):
  Identifies possible glycan compositions from m/z values.
- [`comp_to_struc()`](https://glycoverse.github.io/glyanno/dev/reference/comp_to_struc.md):
  Maps glycan compositions to potential glycan structures.
- [`enhance_comp()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_comp.md):
  Refines generic glycan compositions into concrete,
  monosaccharide-specific compositions.
- [`enhance_struc()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_struc.md):
  Resolves low-resolution glycan structures into high-resolution
  structures with linkage information.

![](diagram.png)

``` r
library(glyanno)
library(glydb)
#> Loading required package: glyrepr
```

**Note:** This package relies on
[glyrepr](https://github.com/glycoverse/glyrepr) and
[glydb](https://github.com/glycoverse/glydb). We recommend familiarizing
yourself with these packages, especially `glyrepr`, before using
`glyanno`. Additionally, a basic understanding of IUPAC-condensed
notation is recommended for this vignette. You can refer to this
[tutorial](https://glycoverse.github.io/glyrepr/articles/iupac.html).

## Workflow: M/z -\> Composition -\> Structure

Let’s demonstrate these functions step-by-step, starting with a single
m/z value.

``` r
mz_to_comp(406.1325, charge = 1, adduct = "Na+")
#> # A tibble: 4 × 2
#>      mz composition    
#>   <dbl> <comp>         
#> 1  406. Gal(1)GalNAc(1)
#> 2  406. Gal(1)GlcNAc(1)
#> 3  406. Glc(1)GlcNAc(1)
#> 4  406. Man(1)GlcNAc(1)
```

This returns every composition in `glydb` matching the m/z of 406.1325.
In practice, returning “all” possibilities can be overwhelming; you will
often want to constrain the search space based on biological context.

`glydb` provides helper functions to filter the database. For example,
you might want to limit the search to only human O-GalNAc glycans.

``` r
my_db <- glydb_compositions(species = "Homo sapiens", glycan_type = "O-GalNAc")
my_db
#> <glycan_composition[124]>
#> [1] Gal(1)GalNAc(1)
#> [2] Gal(1)GlcNAc(1)GalNAc(1)
#> [3] GlcNAc(1)GalNAc(1)
#> [4] GlcNAc(2)GalNAc(1)
#> [5] GalNAc(2)
#> [6] Gal(2)GlcNAc(1)GalNAc(1)Neu5Ac(1)
#> [7] Gal(1)GalNAc(2)Neu5Ac(2)
#> [8] Gal(1)GalNAc(1)Fuc(1)
#> [9] Gal(2)GlcNAc(2)GalNAc(2)Fuc(1)
#> [10] Gal(2)GlcNAc(1)GalNAc(1)
#> ... (114 more not shown)
```

`my_db` is simply a
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
vector. You can pass it to the `db` argument of
[`mz_to_comp()`](https://glycoverse.github.io/glyanno/dev/reference/mz_to_comp.md)
to filter results.

``` r
mz_to_comp(406.1325, charge = 1, adduct = "Na+", db = my_db)
#> # A tibble: 1 × 2
#>      mz composition    
#>   <dbl> <comp>         
#> 1  406. Gal(1)GalNAc(1)
```

The results are now significantly more relevant to our specific context.

With the composition identified, we can proceed to determine potential
structures. The m/z value 406.1325 corresponds to the composition
Gal(1)GalNAc(1). What are the possible structures for this composition?

``` r
struc_db <- glydb_structures(species = "Homo sapiens", glycan_type = "O-GalNAc")
comp_to_struc("Gal(1)GalNAc(1)", db = struc_db)
#> # A tibble: 2 × 2
#>   composition     structure          
#>   <comp>          <struct>           
#> 1 Gal(1)GalNAc(1) Gal(b1-3)GalNAc(a1-
#> 2 Gal(1)GalNAc(1) Gal(a1-3)GalNAc(a1-
```

Note that here we use
[`glydb_structures()`](https://glycoverse.github.io/glydb/reference/glydb_structures.html)
instead of
[`glydb_compositions()`](https://glycoverse.github.io/glydb/reference/glydb_compositions.html)
to create a structure-specific database.

Two structures are possible, one is Core 1, the other is Core 5.
Sometimes we just want a “most possible” result. In this case, you can
set `return_best` to `TRUE`:

``` r
comp_to_struc("Gal(1)GalNAc(1)", db = struc_db, return_best = TRUE)
#> # A tibble: 1 × 2
#>   composition     structure          
#>   <comp>          <struct>           
#> 1 Gal(1)GalNAc(1) Gal(b1-3)GalNAc(a1-
```

Core 1 is more common than Core 5, so it was kept. The function keeps
the structure with most citations for multiple matches.

Note that all functions in `glyanno` works vectorizedly:

``` r
comp_to_struc(c("Gal(1)GalNAc(1)", "GlcNAc(1)GalNAc(1)"), db = struc_db, return_best = TRUE)
#> # A tibble: 2 × 2
#>   composition        structure             
#>   <comp>             <struct>              
#> 1 Gal(1)GalNAc(1)    Gal(b1-3)GalNAc(a1-   
#> 2 GlcNAc(1)GalNAc(1) GlcNAc(b1-3)GalNAc(a1-
```

If you set `return_best` to `TRUE`, a common pattern is to directly
fetch the second column from the result:

``` r
comp_to_struc(c("Gal(1)GalNAc(1)", "GlcNAc(1)GalNAc(1)"), db = struc_db, return_best = TRUE)[[2]]
#> <glycan_structure[2]>
#> [1] Gal(b1-3)GalNAc(a1-
#> [2] GlcNAc(b1-3)GalNAc(a1-
#> # Unique structures: 2
```

This vector always has the same length as the input.

## Enhancing Compositions and Structures

Researchers often need to refine generic annotations into more specific
ones. For instance, you may wish to convert generic compositions into
specific monosaccharide lists, or assign potential linkages to a
topology-only structure.
[`enhance_comp()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_comp.md)
and
[`enhance_struc()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_struc.md)
are designed for this purpose.

Both functions return a tibble with two columns: `raw` (the original
input) and `enhanced` (the potential high-resolution candidates).

``` r
enhance_comp("Hex(1)HexNAc(1)")
#> # A tibble: 4 × 2
#>   raw             enhanced       
#>   <comp>          <comp>         
#> 1 Hex(1)HexNAc(1) Gal(1)GalNAc(1)
#> 2 Hex(1)HexNAc(1) Gal(1)GlcNAc(1)
#> 3 Hex(1)HexNAc(1) Glc(1)GlcNAc(1)
#> 4 Hex(1)HexNAc(1) Man(1)GlcNAc(1)
```

``` r
enhance_struc("Gal(??-?)GalNAc(??-")
#> # A tibble: 9 × 2
#>   raw                 enhanced           
#>   <struct>            <struct>           
#> 1 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(a1-
#> 2 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(b1-
#> 3 Gal(??-?)GalNAc(??- Gal(a1-3)GalNAc(b1-
#> 4 Gal(??-?)GalNAc(??- Gal(b1-4)GalNAc(b1-
#> 5 Gal(??-?)GalNAc(??- Gal(a1-6)GalNAc(a1-
#> 6 Gal(??-?)GalNAc(??- Gal(b1-6)GalNAc(a1-
#> 7 Gal(??-?)GalNAc(??- Gal(b1-6)GalNAc(b1-
#> 8 Gal(??-?)GalNAc(??- Gal(a1-3)GalNAc(a1-
#> 9 Gal(??-?)GalNAc(??- Gal(b1-4)GalNAc(a1-
```

Similarly, providing a custom database will narrow down the results to
biologically relevant candidates.

``` r
enhance_comp("Hex(1)HexNAc(1)", db = glydb_compositions(species = "Homo sapiens", glycan_type = "O-GalNAc"))
#> # A tibble: 1 × 2
#>   raw             enhanced       
#>   <comp>          <comp>         
#> 1 Hex(1)HexNAc(1) Gal(1)GalNAc(1)
```

``` r
enhance_struc("Gal(??-?)GalNAc(??-", db = glydb_structures(species = "Homo sapiens", glycan_type = "O-GalNAc"))
#> # A tibble: 2 × 2
#>   raw                 enhanced           
#>   <struct>            <struct>           
#> 1 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(a1-
#> 2 Gal(??-?)GalNAc(??- Gal(a1-3)GalNAc(a1-
```

You can set `return_best` to `TRUE` as well.

``` r
enhance_struc(
  "Gal(??-?)GalNAc(??-",
  db = glydb_structures(species = "Homo sapiens", glycan_type = "O-GalNAc"),
  return_best = TRUE
)
#> # A tibble: 1 × 2
#>   raw                 enhanced           
#>   <struct>            <struct>           
#> 1 Gal(??-?)GalNAc(??- Gal(b1-3)GalNAc(a1-
```

## Bidirectional Conversion

While `glyanno` focuses on inferring high-resolution information from
low-resolution data, the `glycoverse` ecosystem also provides
complementary functions to perform the reverse operations—simplifying
detailed structures back to their lower-resolution representations.

A summary of these relationships:

| Enhance                                                                                  | Back                                                                                                              |
|------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| [`mz_to_comp()`](https://glycoverse.github.io/glyanno/dev/reference/mz_to_comp.md)       | [`calculate_mz()`](https://glycoverse.github.io/glyanno/dev/reference/calculate_mz.md)                            |
| [`comp_to_struc()`](https://glycoverse.github.io/glyanno/dev/reference/comp_to_struc.md) | [`glyrepr::as_glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/as_glycan_composition.html)   |
| [`enhance_comp()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_comp.md)   | [`glyrepr::convert_to_generic()`](https://glycoverse.github.io/glyrepr/reference/convert_to_generic.html)         |
| [`enhance_struc()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_struc.md) | [`glyrepr::reduce_structure_level()`](https://glycoverse.github.io/glyrepr/reference/reduce_structure_level.html) |
