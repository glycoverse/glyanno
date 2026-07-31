# glyanno (development version)

* `mz_to_comp()`, `enhance_comp()`, `comp_to_struc()`, and `enhance_struc()`
  now include a `confidence` column in expanded results for customized
  prioritization. (#23)

# glyanno 0.6.0

## New features

* New `enhance_struc_denovo()` reconstructs topological N-glycans and falls back to the topological database when reconstruction is not possible. (#16, #17, #18)

## Breaking changes

* Remove `fill_anomer_pos()` from the exported namespace. Use `glyrepr::fill_anomer_pos()` instead. (#15)

## Minor improvements and bug fixes

* `enhance_struc()` now accepts duplicated input structures, preserves their order and per-occurrence results, and matches each unique non-missing structure only once. (#21)
* `calculate_mz()`, `comp_to_struc()`, `enhance_comp()`, and `mz_to_comp()` now process vector inputs substantially faster by reusing prepared databases and direct lookups. (#20)
* `struc_to_glytoucan()` now looks up accessions from `glydb::glydb_data` first and only falls back to the online GlycanFormatConverter API for missing structures. (#14)

# glyanno 0.5.0

## New features

* Add `struc_to_glytoucan()` to map glycan structures to GlyTouCan accessions (#11).
* Add `fill_anomer_pos()` to fill in missing anomeric position information in glycan structures (#13).

# glyanno 0.4.2

## Minor improvements and bug fixes

* Fix a bug in `enhance_struc()` introduced by the breaking change in glyrepr 0.11.0 (#9).
* Update vignette to mention that you can use a structure database with different structure levels.

# glyanno 0.4.1

## Minor improvements and bug fixes

* `comp_to_struc()` now explicitly returns when the output is empty.

# glyanno 0.4.0

## Breaking changes

* `return_best = TRUE` now returns a vector with the same length as input (with `NA` for unmatched glycans) instead of a tibble with potentially fewer rows. Affected functions: `enhance_comp()`, `comp_to_struc()`, `mz_to_comp()`, `enhance_struc()`.

## Minor improvements and bug fixes

* `comp_to_struc()` now errors when concrete compositions are matched against a generic structure database, instead of silently returning empty results.
* `db` paramter cannot be of zero length now, in which case an error will be raised.

# glyanno 0.3.0

## Breaking changes

* Remove `to_level` parameter from `enhance_struc()`.

# glyanno 0.2.0

## New features

* Add `return_best` parameter to keep the most cited match when multiple compositions or structures are possible.

# glyanno 0.1.2

## Minor improvements and bug fixes

* Update dependency strategy to use the r-universe repo.

# glyanno 0.1.1

## Minor improvements and bug fixes

* Adapt to glyrepr 0.10.0.

# glyanno 0.1.0

* First GitHub release.
