# Changelog

## glyanno (development version)

### New features

- [`enhance_struc()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_struc.md)
  gains `method = "denovo"` to reconstruct topological complex and
  hybrid N-glycans from branch candidates and high-mannose N-glycans
  from biosynthetic precursor subtrees. Database fallbacks are
  normalized to topological level. (#16)

### Breaking changes

- Remove
  [`fill_anomer_pos()`](https://glycoverse.github.io/glyrepr/reference/fill_anomer_pos.html)
  from the exported namespace. Use
  [`glyrepr::fill_anomer_pos()`](https://glycoverse.github.io/glyrepr/reference/fill_anomer_pos.html)
  instead. (#15)

### Minor improvements and bug fixes

- [`struc_to_glytoucan()`](https://glycoverse.github.io/glyanno/dev/reference/struc_to_glytoucan.md)
  now looks up accessions from
  [`glydb::glydb_data`](https://glycoverse.github.io/glydb/reference/glydb_data.html)
  first and only falls back to the online GlycanFormatConverter API for
  missing structures. (#14)

## glyanno 0.5.0

### New features

- Add
  [`struc_to_glytoucan()`](https://glycoverse.github.io/glyanno/dev/reference/struc_to_glytoucan.md)
  to map glycan structures to GlyTouCan accessions (#11).
- Add
  [`fill_anomer_pos()`](https://glycoverse.github.io/glyrepr/reference/fill_anomer_pos.html)
  to fill in missing anomeric position information in glycan structures
  (#13).

## glyanno 0.4.2

### Minor improvements and bug fixes

- Fix a bug in
  [`enhance_struc()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_struc.md)
  introduced by the breaking change in glyrepr 0.11.0 (#9).
- Update vignette to mention that you can use a structure database with
  different structure levels.

## glyanno 0.4.1

### Minor improvements and bug fixes

- [`comp_to_struc()`](https://glycoverse.github.io/glyanno/dev/reference/comp_to_struc.md)
  now explicitly returns when the output is empty.

## glyanno 0.4.0

### Breaking changes

- `return_best = TRUE` now returns a vector with the same length as
  input (with `NA` for unmatched glycans) instead of a tibble with
  potentially fewer rows. Affected functions:
  [`enhance_comp()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_comp.md),
  [`comp_to_struc()`](https://glycoverse.github.io/glyanno/dev/reference/comp_to_struc.md),
  [`mz_to_comp()`](https://glycoverse.github.io/glyanno/dev/reference/mz_to_comp.md),
  [`enhance_struc()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_struc.md).

### Minor improvements and bug fixes

- [`comp_to_struc()`](https://glycoverse.github.io/glyanno/dev/reference/comp_to_struc.md)
  now errors when concrete compositions are matched against a generic
  structure database, instead of silently returning empty results.
- `db` paramter cannot be of zero length now, in which case an error
  will be raised.

## glyanno 0.3.0

### Breaking changes

- Remove `to_level` parameter from
  [`enhance_struc()`](https://glycoverse.github.io/glyanno/dev/reference/enhance_struc.md).

## glyanno 0.2.0

### New features

- Add `return_best` parameter to keep the most cited match when multiple
  compositions or structures are possible.

## glyanno 0.1.2

### Minor improvements and bug fixes

- Update dependency strategy to use the r-universe repo.

## glyanno 0.1.1

### Minor improvements and bug fixes

- Adapt to glyrepr 0.10.0.

## glyanno 0.1.0

- First GitHub release.
