# Manage the persistent annotation cache

`build_glyanno_cache()` is the magic that makes `glyanno` fast. It
precalculates database views and annotation indexes from the installed
`glydb`. You will only need to call this function once on your computer.
When `glydb` is updated, you'll also need to call this function to
update the cache. But don't worry. We'll remind you then.

- `build_glyanno_cache()` reuses a valid cache unless `force = TRUE`.

- `glyanno_cache_info()` inspects it without building it, and

- `clear_glyanno_cache()` removes it and clears the session cache.

## Usage

``` r
build_glyanno_cache(force = FALSE)

glyanno_cache_info()

clear_glyanno_cache()
```

## Arguments

- force:

  Rebuild even when the cache is current.

## Value

`build_glyanno_cache()` invisibly returns the cache file path.
`glyanno_cache_info()` returns a list with `path`, `status` (one of
`"missing"`, `"outdated"`, `"unreadable"`, or `"ready"`), `size_bytes`,
`expected`, and `stored` version metadata. `clear_glyanno_cache()`
invisibly returns `NULL`.

## Details

The cache lives in `tools::R_user_dir("glyanno", "cache")`. Set the
`glyanno.cache_dir` option to use another directory. Only an explicit
call to `build_glyanno_cache()` writes persistent data. Without a disk
cache, annotation builds the required views in memory for the current
session.

On interactive attachment, a startup message reminds you to build a
missing, outdated, or unreadable cache. No computation or questions
block attachment. The cache is invalidated by a changed glydb version,
glyrepr or igraph version, R major/minor version, collation locale, or
internal cache schema. Restart R after updating dependencies.
Development data changes without a version bump require `force = TRUE`.

Cached data include concrete/generic compositions, intact/topological
structures, confidence ranks, exact and compatible composition indexes,
residue counts, species/type membership, floating-structure flags,
masses for all built-in dictionaries, and local GlyTouCan accession
keys. Structure vectors retain their parsed graphs. Arbitrary structure
matching still uses glymotif at query time. Custom databases are never
persisted.

A complete replacement is written before the previous cache is replaced;
failed builds preserve the previous file. Only one completed cache is
kept. Metadata inspection reads a small file header and checks the file
size. The compressed annotation payload is loaded and validated on first
use. Use `force = TRUE` to repair a damaged payload whose header is
still valid.

## Examples

``` r
glyanno_cache_info()
#> $path
#> [1] "/home/runner/.cache/R/glyanno/annotation-v1.rds"
#> 
#> $status
#> [1] "missing"
#> 
#> $size_bytes
#> [1] 0
#> 
#> $expected
#> $expected$schema
#> [1] 2
#> 
#> $expected$glydb
#> [1] "0.7.0"
#> 
#> $expected$glyrepr
#> [1] "1.1.0"
#> 
#> $expected$igraph
#> [1] "2.3.3"
#> 
#> $expected$R
#> [1] "4.6"
#> 
#> $expected$collation
#> [1] "C"
#> 
#> 
#> $stored
#> NULL
#> 
```
